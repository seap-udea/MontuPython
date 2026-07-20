###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import ephem as pyephem
import pandas as pd
import numpy as np

from tabulate import tabulate
from functools import lru_cache
from scipy.optimize import brentq

# Planet classes
from pymeeus.Mercury import Mercury as pymeeus_Mercury
from pymeeus.Venus import Venus as pymeeus_Venus
from pymeeus.Mars import Mars as pymeeus_Mars
from pymeeus.Jupiter import Jupiter as pymeeus_Jupiter
from pymeeus.Saturn import Saturn as pymeeus_Saturn
from pymeeus.Uranus import Uranus as pymeeus_Uranus
from pymeeus.Neptune import Neptune as pymeeus_Neptune
from pymeeus.Moon import Moon as pymeeus_Moon
from pymeeus.Sun import Sun as pymeeus_Sun

###############################################################
# Module constants
###############################################################
# Planets
PLANETARY_IDS = dict(
    SUN = 10,
    MERCURY = 1,
    VENUS = 2,
    EARTH = 399,
    MOON = 301,
    MARS = 4,
    JUPITER = 5,
    SATURN = 6,
    URANUS = 7,
    NEPTUNE = 8,
)
PLANETARY_NAMES = {str(v): k for k, v in PLANETARY_IDS.items()}

# Name of quarters of the moon in pymeeus
PYMEEUS_QUARTERS = ['new','first','full','last']

# Field metadata for show_position / show_conditions (key, label, unit)
_SKY_POSITION_SPECS = (
    ('Name', 'Name', ''),
    ('RAJ2000', 'RA (J2000)', 'h'),
    ('DecJ2000', 'Dec (J2000)', 'deg'),
    ('RAEpoch', 'RA (epoch)', 'h'),
    ('DecEpoch', 'Dec (epoch)', 'deg'),
    ('RAGeo', 'RA (geocentric)', 'h'),
    ('DecGeo', 'Dec (geocentric)', 'deg'),
    ('az', 'Azimuth', 'deg'),
    ('el', 'Elevation', 'deg'),
)

_SEBAU_CONDITION_SPECS = (
    ('Name', 'Name', ''),
    ('ha', 'Hour angle', 'h'),
    ('Vmag', 'Visual magnitude', 'mag'),
    ('rise_time', 'Rise time (UTC)', 'jed'),
    ('rise_az', 'Rise azimuth', 'deg'),
    ('set_time', 'Set time (UTC)', 'jed'),
    ('set_az', 'Set azimuth', 'deg'),
    ('transit_time', 'Transit time (UTC)', 'jed'),
    ('transit_el', 'Transit elevation', 'deg'),
    ('elongation', 'Elongation from Sun', 'deg'),
    ('earth_distance', 'Distance from Earth', 'au'),
    ('sun_distance', 'Distance from Sun', 'au'),
    ('angsize', 'Angular diameter', 'arcsec'),
    ('phase', 'Illuminated fraction', 'percent'),
    ('hlat', 'Heliocentric latitude', 'deg'),
    ('hlon', 'Heliocentric longitude', 'deg'),
    ('hlong', 'Heliocentric longitude (alt.)', 'deg'),
    ('is_circumpolar', 'Circumpolar', 'bool'),
    ('is_neverup', 'Never rises', 'bool'),
)

_STAR_CONDITION_SPECS = (
    ('Name', 'Name', ''),
    ('ha', 'Hour angle', 'h'),
    ('Vmag', 'Visual magnitude', 'mag'),
    ('rise_time', 'Rise time (UTC)', 'jed'),
    ('rise_az', 'Rise azimuth', 'deg'),
    ('set_time', 'Set time (UTC)', 'jed'),
    ('set_az', 'Set azimuth', 'deg'),
    ('transit_time', 'Transit time (UTC)', 'jed'),
    ('transit_el', 'Transit elevation', 'deg'),
    ('elongation', 'Elongation from Sun', 'deg'),
    ('is_circumpolar', 'Circumpolar', 'bool'),
    ('is_neverup', 'Never rises', 'bool'),
)

_STAR_POSITION_SPECS = (
    ('Name', 'Name', ''),
    ('RAJ2000', 'RA (J2000)', 'h'),
    ('DecJ2000', 'Dec (J2000)', 'deg'),
    ('RAJ2000t', 'RA (J2000, proper motion)', 'h'),
    ('DecJ2000t', 'Dec (J2000, proper motion)', 'deg'),
    ('RAEpoch', 'RA (epoch)', 'h'),
    ('DecEpoch', 'Dec (epoch)', 'deg'),
    ('RAGeo', 'RA (geocentric)', 'h'),
    ('DecGeo', 'Dec (geocentric)', 'deg'),
    ('az', 'Azimuth', 'deg'),
    ('el', 'Elevation', 'deg'),
)


def _format_sky_value(key, value, unit):
    """Format one sky-quantity for human-readable reports."""
    if value is None:
        return '—'
    if unit == 'bool':
        return 'yes' if value else 'no'
    if unit == '':
        return str(value)
    if unit == 'h':
        return f"{montu.D2S(float(value))} h"
    if unit == 'deg':
        return f"{float(value):.6f}°"
    if unit == 'mag':
        return f"{float(value):.2f} mag"
    if unit == 'au':
        return f"{float(value):.6f} AU"
    if unit == 'arcsec':
        return f"{float(value):.3f}\""
    if unit == 'percent':
        return f"{float(value):.2f} %"
    if unit == 'jed':
        if float(value) == 0:
            return '—'
        return montu.Time(float(value), format='jd').readable.datemix
    return str(value)

###############################################################
# Class Sebau
###############################################################
class Sebau(object):
    """Base class for a solar-system celestial body (Sun, Moon, planet).

    The name *sebau* (Egyptian: *sbꜣw*, 'star') is used here as a generic
    term for any bright celestial object tracked across the sky.

    This class is not intended to be instantiated directly; use the
    subclasses :class:`Sun`, :class:`Moon`, and :class:`Planet` instead.

    Attributes
    ----------
    position : list or montu.Dictobj
        When ``store=False`` (single call), a :class:`montu.Dictobj` with
        the latest equatorial and horizontal coordinates.  When
        ``store=True`` (accumulating calls), a list of position dicts that
        can be converted to a DataFrame with :meth:`tabulate_positions`.
    condition : list or montu.Dictobj
        Same accumulation pattern as *position* but containing
        observational conditions (rise/set/transit times, elongation, etc.).
    """

    def __init__(self):
        # Basic attributes
        self.reset_store()

    @staticmethod
    def _observer_copy(observer, date=None):
        """Create a fresh PyEphem observer with the same site properties."""
        site = pyephem.Observer()
        site.lon = observer.site.lon
        site.lat = observer.site.lat
        site.pressure = observer.site.pressure
        site.temp = observer.site.temp
        site.elevation = observer.site.elevation
        if date is not None:
            site.date = date
        else:
            site.date = observer.site.date
        return site

    def _remember_sky_context(self, at, observer):
        """Store the epoch and site used in the latest sky computation."""
        if at is None:
            at = montu.Time()
        self._sky_at = at
        self._sky_observer = observer

    def _sky_body_label(self):
        return getattr(self, 'name', None) or getattr(self.seba, 'name', 'Object')

    def _sky_position_specs(self):
        return _SKY_POSITION_SPECS

    def _sky_condition_specs(self):
        return _SEBAU_CONDITION_SPECS

    def _get_latest_sky_record(self, attr):
        """Return the latest position or condition dict, if any."""
        value = getattr(self, attr, None)
        if isinstance(value, montu.Dictobj):
            return dict(value.__dict__)
        if isinstance(value, list) and value:
            row = value[-1]
            return dict(row) if isinstance(row, dict) else row
        if isinstance(value, pd.DataFrame) and not value.empty:
            return value.iloc[-1].to_dict()
        return None

    def _format_sky_header(self, *, record):
        """Build epoch and site lines shared by position/condition reports."""
        at = getattr(self, '_sky_at', None)
        if at is None and record and record.get('jed') is not None:
            at = montu.Time(float(record['jed']), format='jd')
        observer = getattr(self, '_sky_observer', None)

        lines = []
        if at is not None:
            at.readable  # populate calendar strings on demand
            lines.append(
                f"  Epoch: {at.readable.datemix} / {at.readable.datepro}  "
                f"(JED {at.jed:.6f})"
            )
        else:
            lines.append("  Epoch: (unknown)")

        if observer is not None:
            site_label = observer.site_name or "Custom site"
            if observer.site_id:
                site_label += f" [{observer.site_id}]"
            lines.append(
                f"  Site: {site_label} — lat {observer.lat:.6f}°, "
                f"lon {observer.lon:.6f}°, {observer.height * 1000:.0f} m  "
                f"(P={observer.pressure} mbar, T={observer.temperature} °C)"
            )
        else:
            lines.append("  Site: (unknown)")
        return lines

    def _format_sky_record(self, record, specs):
        lines = []
        for key, label, unit in specs:
            if key not in record:
                continue
            value = _format_sky_value(key, record.get(key), unit)
            lines.append(f"  {label}: {value}")
        return lines

    def show_position(self):
        """Print the latest :meth:`where_in_sky` result with epoch, site, and units.

        If :meth:`where_in_sky` has not been called yet, prints a short notice
        instead.

        Examples
        --------
        >>> import montu
        >>> spica = montu.Star(montu.Stars(subset='bright', ProperName='Spica').data.iloc[0])
        >>> thebes = montu.Observer(site='thebes')
        >>> mtime = montu.Time('bce 1500-06-21 12:00:00', calendar='mixed')
        >>> spica.where_in_sky(at=mtime, observer=thebes)
        >>> spica.show_position()  # doctest: +SKIP
        """
        record = self._get_latest_sky_record('position')
        label = self._sky_body_label()
        if record is None:
            message = (
                f"{label} — no sky position stored.\n"
                "  Call where_in_sky(at=..., observer=...) first."
            )
            print(message)
            return

        lines = [f"{label} — sky position"]
        lines.extend(self._format_sky_header(record=record))
        lines.extend(self._format_sky_record(record, self._sky_position_specs()))
        print("\n".join(lines))

    def show_conditions(self):
        """Print the latest :meth:`conditions_in_sky` result with epoch, site, and units.

        Planet, Sun, and Moon reports include heliocentric distances and phase;
        :class:`montu.Star` reports omit planet-specific fields.

        If :meth:`conditions_in_sky` has not been called yet, prints a short
        notice instead.

        Examples
        --------
        >>> import montu
        >>> mars = montu.Planet('Mars')
        >>> giza = montu.Observer(lon=31.134, lat=29.979, height=0.075)
        >>> mtime = montu.Time('-1000-03-21 20:00:00')
        >>> mars.conditions_in_sky(at=mtime, observer=giza)
        >>> mars.show_conditions()  # doctest: +SKIP
        """
        record = self._get_latest_sky_record('condition')
        label = self._sky_body_label()
        if record is None:
            message = (
                f"{label} — no sky conditions stored.\n"
                "  Call conditions_in_sky(at=..., observer=...) first."
            )
            print(message)
            return

        lines = [f"{label} — sky conditions"]
        lines.extend(self._format_sky_header(record=record))
        lines.extend(self._format_sky_record(record, self._sky_condition_specs()))
        print("\n".join(lines))

    def _observer_events(self, observer):
        """Compute rise/set/transit event data without deprecated body attrs."""
        event_site = self._observer_copy(observer)
        is_circumpolar = False
        is_neverup = False

        rise_time = 0
        rise_az = 2 * np.pi
        set_time = 0
        set_az = 2 * np.pi
        transit_time = 0
        transit_alt = 2 * np.pi

        try:
            rise_date = event_site.next_rising(self.seba)
            rise_time = float(rise_date)
            rise_site = self._observer_copy(observer, rise_date)
            self.seba.compute(rise_site)
            rise_az = float(self.seba.az)
        except pyephem.AlwaysUpError:
            is_circumpolar = True
        except pyephem.NeverUpError:
            is_neverup = True

        try:
            set_date = event_site.next_setting(self.seba)
            set_time = float(set_date)
            set_site = self._observer_copy(observer, set_date)
            self.seba.compute(set_site)
            set_az = float(self.seba.az)
        except pyephem.AlwaysUpError:
            is_circumpolar = True
        except pyephem.NeverUpError:
            is_neverup = True

        try:
            transit_date = event_site.next_transit(self.seba)
            transit_time = float(transit_date)
            transit_site = self._observer_copy(observer, transit_date)
            self.seba.compute(transit_site)
            transit_alt = float(self.seba.alt)
        except (pyephem.AlwaysUpError, pyephem.NeverUpError):
            pass

        # Restore body to the original observer epoch after temporary event computations.
        self.seba.compute(observer.site)

        return dict(
            rise_time=rise_time,
            rise_az=rise_az,
            set_time=set_time,
            set_az=set_az,
            transit_time=transit_time,
            transit_alt=transit_alt,
            is_circumpolar=is_circumpolar,
            is_neverup=is_neverup,
        )

    def where_in_sky(self,at=None,observer=None,store=False):
        """Compute the current position of the body in the sky.

        Updates ``self.position`` with equatorial (J2000 and epoch) and
        horizontal coordinates at the requested time and place.

        Parameters
        ----------
        at : montu.Time
            Epoch of the observation.
        observer : montu.Observer
            Observing site.
        store : bool, optional
            If ``False`` (default), ``self.position`` is replaced with a
            :class:`montu.Dictobj` containing the latest result.
            If ``True``, the result is appended to ``self.position`` (a list)
            for later bulk conversion with :meth:`tabulate_positions`.

        Examples
        --------
        >>> import montu
        >>> giza = montu.Observer(lon=31.134, lat=29.979, height=0.075)
        >>> mtime = montu.Time('-1000-03-21 06:00:00')
        >>> sun = montu.Sun()
        >>> sun.where_in_sky(at=mtime, observer=giza)
        >>> print(sun.position.az, sun.position.el)
        """
        if at is None:
            at = montu.Time()
        self._remember_sky_context(at, observer)
        self._compute_ephemerides(at.jed,observer)

        # Basic store
        position = {
            'tt':int(at.tt),'jed':at.jed,
            'Name':self.seba.name,'RAJ2000':self.seba.a_ra*montu.RAD/15,
            'DecJ2000':self.seba.a_dec*montu.RAD,'RAEpoch':self.seba.ra*montu.RAD/15,
            'DecEpoch':self.seba.dec*montu.RAD,'RAGeo':self.seba.g_ra*montu.RAD/15,
            'DecGeo':self.seba.g_dec*montu.RAD,'el':self.seba.alt*montu.RAD,
            'az':self.seba.az*montu.RAD,
        }
        # Accumulating store
        self.position_store = store
        if store:
            self.position += [position]
        else:
            self.position = montu.Dictobj(dict=position)
    
    def conditions_in_sky(self,at=None,observer=None,store=False):
        """Compute full observational conditions for the body.

        Extends :meth:`where_in_sky` with rise/set/transit times, elongation,
        angular size, phase, heliocentric coordinates, and more.

        Parameters
        ----------
        at : montu.Time
            Epoch of the observation.
        observer : montu.Observer
            Observing site.
        store : bool, optional
            If ``False`` (default), replace ``self.condition`` with a fresh
            :class:`montu.Dictobj`.  If ``True``, append to ``self.condition``
            list for later batch conversion.

        Examples
        --------
        >>> import montu
        >>> giza = montu.Observer(lon=31.134, lat=29.979, height=0.075)
        >>> mtime = montu.Time('-1000-03-21 18:00:00')
        >>> mars = montu.Planet('Mars')
        >>> mars.conditions_in_sky(at=mtime, observer=giza)
        >>> print(mars.condition.rise_time, mars.condition.set_time)
        """
        # First compute
        self.where_in_sky(at,observer,store)
        events = self._observer_events(observer)
        
        # Store
        condition = {
            'tt':int(at.tt),'jed':at.jed,
            'Name':self.seba.name,
            'ha':self.seba.ha*montu.RAD/15,
            'Vmag':self.seba.mag,
            'rise_time':events['rise_time'] + montu.PYEPHEM_JD_REF,
            'rise_az':events['rise_az']*montu.RAD,
            'set_time':events['set_time'] + montu.PYEPHEM_JD_REF,
            'set_az':events['set_az']*montu.RAD,
            'transit_time':events['transit_time'] + montu.PYEPHEM_JD_REF,
            'transit_el':events['transit_alt']*montu.RAD,
            'elongation':self.seba.elong*montu.RAD,'earth_distance':self.seba.earth_distance,
            'sun_distance':self.seba.sun_distance,'is_circumpolar':events['is_circumpolar'],
            'is_neverup':events['is_neverup'],'angsize':self.seba.size,
            'phase':self.seba.phase,'hlat':self.seba.hlat*montu.RAD,'hlon':self.seba.hlon*montu.RAD,
            'hlong':self.seba.hlong*montu.RAD,
        }
        
        self.condition_store = store
        if store:
            self.condition += [condition]
        else:
            self.condition = montu.Dictobj(dict=condition)

    def _compute_ephemerides(self,jed=None,observer=None):
        """Call PyEphem to compute ephemerides for the body.

        Parameters
        ----------
        jed : float, optional
            Julian Ephemeris Day (UTC scale). Defaults to current time.
        observer : montu.Observer
            Observing site used as the reference frame.

        Raises
        ------
        ValueError
            If *observer* is not a :class:`montu.Observer` instance.
        """
        if not isinstance(observer,montu.Observer):
            raise ValueError("You must provide a valid montu.Observer")

        if jed is None:
            jed = montu.Time().jed

        # Update observer site
        observer.site.date = jed - montu.PYEPHEM_JD_REF

        # Compute ephemerides
        self.seba.compute(observer.site)

    def when_it_appears(self,at=None,observer=None):
        """Compute the time when the body first appears in the night sky.

        .. note::
            Not yet implemented.
        """
        pass

    def __str__(self):

        # Positions
        output = f"Object {self.seba.name} positions:\n"
        if isinstance(self.position,montu.Dictobj):
            tabulable = [self.position.__dict__]
        else:
            tabulable = self.position
        output += f"{tabulate(tabulable,headers='keys',tablefmt='github')}'"

        # Conditions
        output += f"\nObject {self.seba.name} conditions:\n"
        if isinstance(self.condition,montu.Dictobj):
            tabulable = [self.condition.__dict__]
        else:
            tabulable = self.condition
        output += f"{tabulate(tabulable,headers='keys',tablefmt='github')}'"
        
        return output
    
    def __repr__(self):
        return self.__str__()

    def reset_store(self):
        """Reset accumulated position and condition lists to empty.

        Call this before starting a new set of accumulating
        ``store=True`` calls.

        Examples
        --------
        >>> sun = montu.Sun()
        >>> sun.reset_store()
        """
        self.position = []
        self.position_store = False
        self.condition = []
        self.condition_store = False
        self._sky_at = None
        self._sky_observer = None

    def tabulate_store(self):
        """Convert accumulated position/condition lists to DataFrames in place.

        Converts ``self.position`` and/or ``self.condition`` from lists of
        dicts into :class:`pandas.DataFrame` objects.
        Call this after finishing a series of ``store=True`` computations.

        Examples
        --------
        >>> import montu
        >>> giza = montu.Observer(lon=31.134, lat=29.979)
        >>> sun = montu.Sun()
        >>> times = [montu.Time('-1000-01-01') + i * montu.DAY for i in range(30)]
        >>> for t in times:
        ...     sun.where_in_sky(at=t, observer=giza, store=True)
        >>> sun.tabulate_store()
        >>> sun.position.head()  # now a DataFrame
        """
        if self.position_store:
            self.position = pd.DataFrame(self.position)
        if self.condition_store:
            self.condition = pd.DataFrame(self.condition)

    def tabulate_positions(self):
        """Convert the position list to a :class:`pandas.DataFrame`.

        Returns
        -------
        pandas.DataFrame
            DataFrame with one row per :meth:`where_in_sky` call made with
            ``store=True``.
        """
        return pd.DataFrame(self.position)
    
    def tabulate_conditions(self):
        """Convert the condition list to a :class:`pandas.DataFrame`.

        Returns
        -------
        pandas.DataFrame
            DataFrame with one row per :meth:`conditions_in_sky` call made
            with ``store=True``.
        """
        return pd.DataFrame(self.condition)
    
    def tabulate_ephemerides(self):
        """Merge position and condition DataFrames into one ephemerides table.

        Calls :meth:`tabulate_store` first; then merges ``self.position``
        and ``self.condition`` on ``['tt', 'jed', 'Name']`` and stores the
        result in ``self.ephemerides``.

        Examples
        --------
        >>> sun.tabulate_ephemerides()
        >>> sun.ephemerides.head()
        """
        self.tabulate_store()

        # Create a unique ephemerides pandas object
        self.ephemerides = self.position
        if self.condition_store:
            self.ephemerides = pd.merge(self.ephemerides,
                                        self.condition,
                                        on=['tt','jed','Name'])
            
    @staticmethod
    def where_in_sky_all_planets(at=None,observer=None):
        pass
        
###############################################################
# Sun Class
###############################################################
class Sun(Sebau):
    """The Sun as a celestial body.

    Inherits all methods from :class:`Sebau`. Provides additional static
    methods for computing solstices, equinoxes, and twilight times.

    Examples
    --------
    >>> import montu
    >>> giza = montu.Observer(lon=31.134, lat=29.979, height=0.075)
    >>> mtime = montu.Time('-1000-03-21 06:00:00')
    >>> sun = montu.Sun()
    >>> sun.where_in_sky(at=mtime, observer=giza)
    >>> print(f"Az={sun.position.az:.2f} deg, El={sun.position.el:.2f} deg")
    """
    def __init__(self):
        super().__init__()
        self.seba = pyephem.Sun()
        self.name = self.seba.name

    def where_in_sky(self, at=None, observer=None, store=False):
        super().where_in_sky(at, observer, store)

    def conditions_in_sky(self, at=None, observer=None, store=False):
        super().conditions_in_sky(at, observer, store)
        
    @staticmethod
    def next_seasons(at=None):
        """Return the Julian days of the four upcoming astronomical seasons.

        Parameters
        ----------
        at : montu.Time
            Reference date from which to search forward.

        Returns
        -------
        vernal_jed : float
            Julian day of the next vernal equinox.
        summer_jed : float
            Julian day of the next summer solstice.
        autumnal_jed : float
            Julian day of the next autumnal equinox.
        winter_jed : float
            Julian day of the next winter solstice.

        Examples
        --------
        >>> import montu
        >>> ve, ss, ae, ws = montu.Sun.next_seasons(at=montu.Time('-1000-01-01'))
        >>> montu.Time(ve, format='jd').get_readable().readable.datepro
        '-1000-03-...'
        """
        date = pyephem.Date(at.jed - montu.PYEPHEM_JD_REF)
        vernal_jed = pyephem.next_vernal_equinox(date) + montu.PYEPHEM_JD_REF
        summer_jed = pyephem.next_summer_solstice(date) + montu.PYEPHEM_JD_REF
        auttumnal_jed = pyephem.next_autumnal_equinox(date) + montu.PYEPHEM_JD_REF
        winter_jed = pyephem.next_winter_solstice(date) + montu.PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed
    
    @staticmethod
    def previous_seasons(at=None):
        """Return the Julian days of the four most recent astronomical seasons.

        Parameters
        ----------
        at : montu.Time
            Reference date from which to search backward.

        Returns
        -------
        vernal_jed : float
            Julian day of the preceding vernal equinox.
        summer_jed : float
            Julian day of the preceding summer solstice.
        autumnal_jed : float
            Julian day of the preceding autumnal equinox.
        winter_jed : float
            Julian day of the preceding winter solstice.

        Examples
        --------
        >>> import montu
        >>> ve, ss, ae, ws = montu.Sun.previous_seasons(at=montu.Time('-1000-06-01'))
        >>> montu.Time(ve, format='jd').get_readable().readable.datepro
        '-1000-03-...'
        """
        date = pyephem.Date(at.jed - montu.PYEPHEM_JD_REF)
        vernal_jed = pyephem.previous_vernal_equinox(date) + montu.PYEPHEM_JD_REF
        summer_jed = pyephem.previous_summer_solstice(date) + montu.PYEPHEM_JD_REF
        auttumnal_jed = pyephem.previous_autumnal_equinox(date) + montu.PYEPHEM_JD_REF
        winter_jed = pyephem.previous_winter_solstice(date) + montu.PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed

    @staticmethod
    def when_is_twilight(day=None, observer=None, sunbelow=-6):
        """Time of start and end of night time (between twilights).

        Parameters:
            at: montu.Time, default = None:
                Time at which the twilight is calculated.

            observer: montu.Observer, default = None:
                Observer which see the object.

            sunbelow: float [deg], default = -6:
                Angle below the horizon on which the astronomical
                twilight is defined.
                For convention: 
                    sunbelow = -6 for civil twilight
                    sunbelow = -12 for nautical twilight
                    sunbelow = -18 for astronomical twilight

        Return:
            dusk_time, down_time: float [julian day]:
                Time of start and end of night: time of astronomical dusk
                time of astronomical down.

            appearance_function: function(Vmag):
                It gives you a function to compute the time when the 
                object starts to be observed and when it dissapears.
                This time depends on the angle below the horizon from which 
                we define the object will be visible under clear sky conditions.

        References:
            https://en.wikipedia.org/wiki/Twilight
        """
        # Get rise and set time for Sun
        sun = montu.Sun()
        sun.conditions_in_sky(at=day,observer=observer)
        set_time = sun.condition.set_time
        rise_time = sun.condition.rise_time
        
        # Routine for calculating twilight
        def is_sun_elevation_at(dt,ref_time,elevation):
            time = ref_time + dt/montu.DAY
            sun._compute_ephemerides(time,observer=observer)
            sun_elevation = sun.seba.alt*montu.RAD
            return sun_elevation - elevation

        # Calculate dusk and dawn time
        xtol = 1 # Tolerance set to 1 second to reduce computing time 
        dusk_time = rise_time + brentq(is_sun_elevation_at,-6*montu.HOUR,0,args=(rise_time,sunbelow),xtol=xtol)/montu.DAY
        dawn_time = set_time + brentq(is_sun_elevation_at,0,+6*montu.HOUR,args=(set_time,sunbelow),xtol=xtol)/montu.DAY

        return dusk_time,dawn_time

###############################################################
# Moon Class
###############################################################
class Moon(Sebau):
    """The Moon as a celestial body.

    Inherits all methods from :class:`Sebau`. Provides additional static
    methods for computing moon phases and quarters.

    Examples
    --------
    >>> import montu
    >>> giza = montu.Observer(lon=31.134, lat=29.979, height=0.075)
    >>> mtime = montu.Time('-1000-03-21 20:00:00')
    >>> moon = montu.Moon()
    >>> moon.conditions_in_sky(at=mtime, observer=giza)
    >>> print(f"Phase: {moon.condition.moon_phase:.2f}")
    """
    def __init__(self):
        super().__init__()
        self.seba = pyephem.Moon()
        self.name = self.seba.name

    def where_in_sky(self, at=None, observer=None, store=False):
        super().where_in_sky(at, observer, store)

    def conditions_in_sky(self, at=None, observer=None, store=False):
        super().conditions_in_sky(at, observer, store)
        
        # Additional fields
        additional_conditions = {
            'colong':[self.seba.colong*montu.RAD],'libration_lat':[self.seba.libration_lat*montu.RAD],
            'libration_long':[self.seba.libration_long*montu.RAD],'moon_phase':[self.seba.moon_phase],
        }

    @staticmethod
    def next_moon_quarters(since=None,
                           starting_at=None,
                           numquarters=4,
                           output='jed',
                           format='dict'):
        """
        Compute the next moon quarters starting at a given date.

        Parameters:
            since: montu.Time, default = None:
                Reference date. If None, mtime is now.

            starting_at: string, default = None:
                Quarter on which the sequence starts. If None the sequence
                will start at the first quarter found afer 'since'.

            numquarters: float, default = 4:
                Number of quarters to find. 

            output: string, default = 'jed':
                Which output do you wante:
                    'jed': julian days.
                    'date': date as a list of values.
                    'mtime': montu.Time object.
                    'datepro': date string in proleptic calendar.
                    'datemix': date string in mixed calendar.

            format: string, default = 'dict':
                which format do you want the return:
                    'dict':
                        dict(full1=[<date>,delta-t,delta-0],
                             last1=[<date>,delta-t,delta-0],
                             ...
                            )
                    'columns':
                        [
                            {'Quarter':'full','Datetime':<date>,'delta-t':<delta-t>,'delta-0':<delta-0>},
                            {'Quarter':'last','Datetime':<date>,'delta-t':<delta-t>,'delta-0':<delta-0>},
                        ]

        Return:
            Dates of quarters.

        Examples:
            Determine the date of quarters of the first synodic month after 2000-01-01
                montu.Moon.next_moon_quarters(since=montu.Time('2000-01-01'),starting_at='new')
        """
        if since is None:
            since = montu.Time()
        
        quarter_dates = dict() if format=='dict' else []
        pyepoch = montu.pymeeus_Epoch(since.jed)

        # Prepare search
        new_pyepoch = montu.pymeeus_Moon.moon_phase(pyepoch-15,target='new')
        next_pyepoch = new_pyepoch
        prev_pyepoch = next_pyepoch

        nquarters = 0
        n = 1
        while nquarters<numquarters:
            # Which is the next wuarter
            quarter = PYMEEUS_QUARTERS[n%4]
            if quarter == 'new':
                new_pyepoch = next_pyepoch
            
            # Search for next quarters
            next_pyepoch = montu.pymeeus_Moon.moon_phase(new_pyepoch,target=quarter)
            delta_quarter = float(next_pyepoch - prev_pyepoch)
            
            # Check if quarter is posterior
            if float(next_pyepoch) > since.jed:
                save_quarter = True

                # If this is the first quarter check that it is the first one
                if nquarters == 0 and (starting_at is not None):
                    if quarter != starting_at:
                        save_quarter = False

                # Save quarter if starting_at condition is fulfilled
                if save_quarter:
                    delta_day = float(next_pyepoch - since.jed)
                    # Extract date
                    if output == 'jed':
                        return_value = float(next_pyepoch)
                    elif output == 'date':
                        return_value = next_pyepoch.get_full_date()
                    elif output == 'mtime':
                        return_value = montu.Time(float(next_pyepoch),format='jd').get_readable()
                    elif output == 'datepro':
                        return_value = montu.Time(float(next_pyepoch),format='jd').get_readable().readable.datepro
                    elif output == 'datemix':
                        return_value = montu.Time(float(next_pyepoch),format='jd').get_readable().readable.datemix
                    else:
                        raise ValueError(f"Output format '{output}' not recognized (valid are 'jed', 'mtime')")
                    
                    # Save quarter
                    if format == 'dict':
                        cycle = int(nquarters/4)+1 if numquarters>4 else ''
                        quarter_dates[quarter+str(cycle)] = [return_value,delta_quarter,delta_day]
                    elif format == 'columns':
                        quarter_dates += [{'Quarter':quarter,'Datetime':return_value,'delta-t':delta_quarter,'delta-0':delta_day}]

                    nquarters += 1
            
            # Sum 
            prev_pyepoch = next_pyepoch
            n += 1
            
        return quarter_dates

def _canonical_planetary_name(name):
    """Return the canonical PLANETARY_IDS key for *name* (case-insensitive)."""
    name_upper = (name or "").strip().upper()
    if name_upper in PLANETARY_IDS:
        return name_upper
    if name_upper in PLANETARY_NAMES:
        return PLANETARY_NAMES[name_upper]
    raise ValueError(
        f"Planet '{name_upper}' not recognized, check variable PLANETARY_NAMES"
    )

###############################################################
# Planet Class
###############################################################
class Planet(Sebau):
    """A solar-system body (Sun, Moon, or a major planet).

    Parameters
    ----------
    name : str
        Body name or NAIF/PyEphem integer ID. Case-insensitive.
        Accepted names: ``'Sun'``, ``'Moon'``, ``'Mercury'``, ``'Venus'``,
        ``'Mars'``, ``'Jupiter'``, ``'Saturn'``, ``'Uranus'``, ``'Neptune'``.

        ``'Sun'`` and ``'Moon'`` return :class:`Sun` and :class:`Moon`
        instances (same objects as :func:`Sun` / :func:`Moon` would build);
        other names return :class:`Planet`.

    Raises
    ------
    ValueError
        If *name* is not a recognised body.

    Attributes
    ----------
    name : str
        Capitalised body name (e.g. ``'Mars'``).
    name_upper : str
        Upper-case body name (e.g. ``'MARS'``).
    id : str
        NAIF body ID as a string (e.g. ``'4'`` for Mars).
    planet_class : type
        Corresponding PyMeeus planet class.

    Examples
    --------
    >>> import montu
    >>> giza = montu.Observer(lon=31.134, lat=29.979, height=0.075)
    >>> mtime = montu.Time('-1000-03-21 20:00:00')
    >>> mars = montu.Planet('Mars')
    >>> mars.where_in_sky(at=mtime, observer=giza)
    >>> print(f"Mars: Az={mars.position.az:.2f}, El={mars.position.el:.2f}")
    >>> isinstance(montu.Planet('Sun'), montu.Sun)
    True
    """
    def __new__(cls, name):
        canonical = _canonical_planetary_name(name)
        if canonical == "SUN":
            return Sun()
        if canonical == "MOON":
            return Moon()
        return super().__new__(cls)

    def __init__(self, name):
        super().__init__()

        canonical = _canonical_planetary_name(name)
        self.name_upper = canonical
        self.id = str(PLANETARY_IDS[canonical])
        self.name = canonical.lower()
        self.name_lower = self.name
        self.name = self.name[0].upper() + self.name[1:]
        self.name_upper = self.name.upper()

        # Find the planet
        exec(f"self.seba = pyephem.{self.name}()")
        self.name = self.seba.name

        # Get class of the planet
        self.planet_class = eval(f"pymeeus_{self.name}")

    def where_in_sky(self, at=None, observer=None, store=False):
        super().where_in_sky(at, observer, store)

    def conditions_in_sky(self, at=None, observer=None, store=False):
        super().conditions_in_sky(at, observer, store)

    def next_planesticies(self,at=None):
        """Compute the two upcoming stations in ecliptic longitude.

        A *planesticio* (Spanish for *stationary point*) is the moment when
        a planet reverses its apparent east--west motion relative to the
        background stars.

        Parameters
        ----------
        at : montu.Time
            Reference epoch from which to search forward.

        Returns
        -------
        jed_station1 : float
            Julian day (UTC) of the first (direct-to-retrograde) station.
        jed_station2 : float
            Julian day (UTC) of the second (retrograde-to-direct) station.

        Examples
        --------
        >>> import montu
        >>> mars = montu.Planet('Mars')
        >>> s1, s2 = mars.next_planesticies(at=montu.Time('-500-01-01'))
        >>> montu.Time(s1, format='jd').get_readable().readable.datepro
        '-500-...-...'
        """
        epoch = montu.pymeeus_Epoch(at.jed)
        
        jed_station1 = self.planet_class.station_longitude_1(epoch)
        jed_station2 = self.planet_class.station_longitude_2(epoch)

        return float(jed_station1),float(jed_station2)