
###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import re
import copy 

import numpy as np
from datetime import datetime
import ephem as pyephem

from scipy.interpolate import interp1d
from pymeeus.Epoch import Epoch as pymeeus_Epoch
from pyplanets.core.epoch import Epoch as pyplanets_Epoch

###############################################################
# Module constants
###############################################################
# Correction for JED
jed_correction_data = np.loadtxt(
    montu.Util._data_path('corrections_dt.dat')
    )
"""
There is an error of between 0.5 and 10 seconds for years 
between 300 c.e. and 1582 c.e. of unknown origin and it seems
to come from the SPICE algorithms or the algorithms calculating
deltat in PyMeeus. This code correct this effect.
"""
JED_CORRECTION = interp1d(
    jed_correction_data[:,0],jed_correction_data[:,1])

# Julian day associated constants
PYEPHEM_JD_REF = 2415020.0
PYEPHEM_MJD_2000 = 36525.0
JED_2000 = 2451545.0
PYMEEUS_JED_2000 = pymeeus_Epoch(JED_2000)

# Astronomical
MIN = 60 # s
HOUR = 60*MIN # s
DAY = 86400 # s
CALYEAR = 365*DAY # s, calendar year
JULYEAR = 365.25*DAY # s, julian year
SIDYEAR = 365.256363004*DAY #s, sidereal year
YEAR = 365.25*DAY # s
CENTURY = 100*YEAR # s
MILLENIUM = 10*YEAR # s

# Abreviation
MONTH_ABREVS = dict(JAN=1,FEB=2,MAR=3,APR=4,MAY=5,JUN=6,
                    JUL=7,AUG=8,SEP=9,OCT=10,NOV=11,DEC=12)

# Day of week (Sunday=1 … Saturday=7); index 0 = Sunday for jed_to_weekday
WEEKDAY_NAMES = ('sunday', 'monday', 'tuesday', 'wednesday',
                 'thursday', 'friday', 'saturday')

# Rounding behavior of time in seconds
"""
Rounding errors may produce strange figures when dealing with 
Julian Days and terrestrial time. We round-up all tt's and 
JD's to avoid this artifacts. If we round to seconds, ie.
ROUND_SECOND_LEVEL = 0, then the rounding level of jd will 
be ROUND_DAY_LEVEL = 1/86400 = 1e-5, i.e. 5.

Rounding routines work in this way: 
    ROUND_SECONDS(1231244.72213) = 1231245.0
"""
ROUND_SECOND_LEVEL = 1 # Decimals included in seconds of date
ROUND_DAY_LEVEL = int(abs(np.log10((ROUND_SECOND_LEVEL+1)/DAY)))+3 
ROUND_SECONDS = lambda seconds:round(seconds,ROUND_SECOND_LEVEL)
ROUND_JULIANDAYS = lambda days:round(days,ROUND_DAY_LEVEL)

# Egyptian civil calendar (sothic)
HORUS_MONTH = dict(I=1,II=2,III=3,IV=4)
HORUS_SEASON = dict(AKHET=1,PERET=2,SHEMU=3,
                    A=1,P=2,S=3,
                    AKH=1,PER=2,SHE=3)
MONTH_HORUS = {str(v): k for k, v in HORUS_MONTH.items()}
SEASON_HORUS = {'0':'Mesut','1':'Akhet','2':'Peret','3':'Shemu'}

# Sothic date: '[hrw HYEAR] MONTH SEASON DAY'  e.g. '[hrw 0] I akhet 1'
SOTHIC_DATE_RE = re.compile(
    r'^\[hrw\s+(-?\d+)\]\s+(\S+)\s+(\S+)\s+(\d+)\s*$',
    re.IGNORECASE,
)
# Mixed-year sothic: '[bce 2026] I akhet 2' or '[-2025] I akhet 2'
SOTHIC_MIXED_DATE_RE = re.compile(
    r'^\[(?P<year_tag>(?!hrw\s)(?:[^\]]+))\]\s+'
    r'(?P<month>\S+)\s+(?P<season>\S+)\s+(?P<day>\d+)\s*$',
    re.IGNORECASE,
)
SOTHIC_TIME_RE = re.compile(
    r'\s(\d{1,2}:\d{2}(?::\d{2})?(?:\.\d+)?)\s*$'
)

# Julian day of the first apokatastasis (coincidence heliakal rise of sopedet and I-Akhet-1)
JED_APOKATASTASIS = 705497.5 # bce 2782-07-20

MJD_EPOCH = np.datetime64('1858-11-17T00:00:00')

def _jed_from_datetime64(dt64):
    """Julian Day (UTC) from a proleptic-Gregorian ``datetime64``."""
    mjd = (dt64 - MJD_EPOCH) / np.timedelta64(1, 'D')
    return float(mjd) + 2400000.5

def _et_from_jed(jed):
    """Ephemeris seconds past J2000 from a Julian Day."""
    return ROUND_SECONDS((jed - JED_2000) * DAY)

def _jed_from_et(et):
    """Julian Day from ephemeris seconds past J2000."""
    return ROUND_JULIANDAYS(JED_2000 + et / DAY)

def _datetime64_from_jed(jed):
    """Proleptic-Gregorian ``datetime64`` from a Julian Day (UTC)."""
    mjd = jed - 2400000.5
    us = int(round(mjd * DAY * 1e6))
    return MJD_EPOCH + np.timedelta64(us, 'us')

def _datestr_from_datetime64(dt64):
    """ISO-like string from a ``datetime64`` timestamp."""
    return str(dt64).replace('T', ' ')


def _days_in_month(year, month, calendar='proleptic', day=1):
    """Days in a month for proleptic-Gregorian or mixed (Julian/Gregorian) rules."""
    if calendar == 'mixed' and pymeeus_Epoch.is_julian(year, month, day):
        if month == 2:
            y = year if year > 0 else year + 1
            return 29 if y % 4 == 0 else 28
        if month in (4, 6, 9, 11):
            return 30
        return 31
    if month in (1, 3, 5, 7, 8, 10, 12):
        return 31
    if month in (4, 6, 9, 11):
        return 30
    y = year
    if y <= 0:
        y += 1
    leap = (y % 4 == 0 and (y % 100 != 0 or y % 400 == 0))
    return 29 if leap else 28


def _apply_civil_offset(y, m, d, h, mi, s,
                        years=0, months=0, days=0, weeks=0,
                        hours=0, minutes=0, seconds=0,
                        calendar='proleptic'):
    """Shift civil date/time components with calendar-aware month lengths."""
    days += weeks * 7
    y = int(y) + years
    m = int(m) + months
    y, m = _normalize_month_year(y, m)
    d = int(d)
    h, mi = int(h), int(mi)
    s = float(s) + seconds

    mi += minutes + int(s // 60)
    s = s % 60
    if s < 0:
        s += 60
        mi -= 1

    h += hours + mi // 60
    mi = mi % 60
    if mi < 0:
        mi += 60
        h -= 1

    d += days + h // 24
    h = h % 24
    if h < 0:
        h += 24
        d -= 1

    while True:
        dim = _days_in_month(y, m, calendar=calendar, day=d)
        if 1 <= d <= dim:
            break
        if d > dim:
            d -= dim
            m += 1
            if m > 12:
                m = 1
                y += 1
        else:
            m -= 1
            if m < 1:
                m = 12
                y -= 1
            d += _days_in_month(y, m, calendar=calendar, day=d)

    sec = int(round(s))
    if sec == 60:
        sec = 0
        mi += 1
        if mi == 60:
            mi = 0
            h += 1
            if h == 24:
                h = 0
                d += 1
    return y, m, d, h, mi, sec


def _mixed_datestr(y, m, d, h, mi, s):
    """Build an ISO-like date string accepted by ``Time(..., calendar='mixed')``."""
    return f'{y}-{int(m):02d}-{int(d):02d} {int(h):02d}:{int(mi):02d}:{int(s):02d}'


def _normalize_month_year(year, month):
    """Carry month overflow into the astronomical year."""
    while month > 12:
        year += 1
        month -= 12
    while month < 1:
        year -= 1
        month += 12
    return year, month


def _sothic_civil_days_to_parts(total_days):
    """Decompose a sothic civil-day offset into years and remainder parts."""
    years, rem = divmod(int(total_days), 365)
    if rem >= 360:
        return years, 0, 0, rem - 360 + 1
    season = rem // 120
    rem -= season * 120
    month = rem // 30
    day = rem % 30 + 1
    return years, month + 1, season + 1, day


class CalendarDelta:
    """Calendar difference between two :class:`Time` instants.

    Returned by :meth:`Time.diff`.  Supports unpacking as
    ``years, days, hours, minutes, seconds`` and conversion helpers
    ``to_days()`` / ``to_years()``.
    """

    __slots__ = ('years', 'months', 'days', 'hours', 'minutes', 'seconds', '_jed_days')

    def __init__(self, years=0, months=0, days=0, hours=0, minutes=0, seconds=0, *, _jed_days=None):
        self.years = int(years)
        self.months = int(months)
        self.days = int(days)
        self.hours = int(hours)
        self.minutes = int(minutes)
        self.seconds = int(seconds)
        self._jed_days = _jed_days

    def __iter__(self):
        return iter((self.years, self.days, self.hours, self.minutes, self.seconds))

    def __neg__(self):
        return CalendarDelta(
            -self.years, -self.months, -self.days,
            -self.hours, -self.minutes, -self.seconds,
            _jed_days=-self._jed_days if self._jed_days is not None else None,
        )

    def __float__(self):
        return self.to_days() * DAY

    def to_days(self):
        """Elapsed Julian days (UTC) between the two instants."""
        if self._jed_days is not None:
            return self._jed_days
        return (
            self.years * CALYEAR + self.months * (CALYEAR / 12)
            + self.days * DAY + self.hours * HOUR + self.minutes * MIN + self.seconds
        ) / DAY

    def to_years(self):
        """Approximate difference in calendar years (365-day civil years)."""
        if self._jed_days is not None:
            return self._jed_days * DAY / CALYEAR
        return self.to_days() * DAY / CALYEAR

    def __repr__(self):
        return (
            f"CalendarDelta(years={self.years}, days={self.days}, "
            f"hours={self.hours}, minutes={self.minutes}, seconds={self.seconds})"
        )


class ReadableTime(montu.Dictobj):
    """Represent all human-readable string representations of a Time object,
    populating them on demand when accessed.
    """
    def __init__(self, time_obj, **kwargs):
        object.__setattr__(self, '_time_obj', time_obj)
        object.__setattr__(self, '_populated', False)
        super().__init__(**kwargs)

    def _ensure_populated(self):
        if not object.__getattribute__(self, '_populated'):
            object.__setattr__(self, '_populated', True)
            object.__getattribute__(self, '_time_obj').get_readable()

    def __getattr__(self, name):
        self._ensure_populated()
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            pass
        raise AttributeError(f"Object of class 'ReadableTime' has no attribute '{name}'")

    def __str__(self):
        self._ensure_populated()
        clean_dict = {k: v for k, v in self.__dict__.items() if not k.startswith('_')}
        return str(clean_dict)

    def __repr__(self):
        return self.__str__()


###############################################################
# Main class
###############################################################
class Time(object):
    """Represent an instant in time, supporting ancient astronomical calendars.

    This is the central class of MontuPython.  It converts dates between the
    proleptic Gregorian calendar, the mixed Julian/Gregorian calendar, the
    ancient Egyptian civil (sothic) calendar, and several numerical
    time-scales used by the underlying ephemeris engines (PyEphem,
    PyMeeus, PyPlanets).

    Parameters
    ----------
    date : str or float
        The date/time to represent.

        When ``format='iso'`` (default), accepted ISO-style strings are::

            '-1000-03-21 06:00:00'          # proleptic BCE: 1001 BCE
            'bce1001-03-21 06:00:00'        # calendar BCE notation
            '1001 b.c.e. 03-21 06:00:00'   # verbose BCE notation
            '2024-05-01 19:00:00'           # CE date

        When ``calendar='sothic'``, the Egyptian civil date format is::

            '[hrw 0] I akhet 1'             # Horus year 0 (≈2782 BCE), I akhet 1
            '[bce 2782] I akhet 1'          # same instant, mixed calendar year (BCE)
            '[-2781] I akhet 1'             # same instant, astronomical mixed year
            '[hrw 1461] I akhet 1 12:00:00' # with civil time of day

        When ``format='jd'``, *date* should be a Julian Day number (float).
        When ``format='tt'``, *date* should be ephemeris seconds past J2000.
    format : {'iso', 'tt', 'jd', 'sothic'}, optional
        Format of *date*.  Default ``'iso'``.
    scale : {'utc', 'tt'}, optional
        Time scale of the input.  ``'utc'`` (default) treats *date* as a
        UTC epoch; ``'tt'`` treats it as Terrestrial Time (TT).
    calendar : {'proleptic', 'mixed', 'sothic'}, optional
        Calendar system for ISO strings. ``'proleptic'`` (default) extends the
        Gregorian calendar before 1582; ``'mixed'`` uses the Julian calendar
        before the Gregorian reform date.
    full : bool, optional
        If ``True``, compute and populate all ``readable`` attributes
        immediately.  Default ``False``.

    Attributes
    ----------
    tt : float
        Ephemeris seconds past J2000 in TT (uniform scale).
    et : float
        Ephemeris seconds past J2000 in UTC.  ``et = tt - deltat``.
    jed : float
        Julian Day in UTC.
    jtd : float
        Julian Day in TT.
    hed : float
        Horus Day in UTC (days since the first apokatastasis, 2782 BCE).
    htd : float
        Horus Day in TT.
    deltat : float
        Difference TT − UTC in seconds.
    bce : bool
        ``True`` if the date is before the common era.
    isjulian : bool
        ``True`` if the date falls within the Julian calendar period.
    readable : montu.Dictobj
        Namespace populated by :meth:`get_readable` with human-readable
        representations:

        * ``datepro`` — proleptic Gregorian string (``'-1000-03-21 06:00:00'``)
        * ``datemix`` — mixed calendar string
        * ``datespice`` — SPICE-format date string
        * ``datesot`` — Egyptian civil (sothic) date string
        * ``comps`` — tuple ``(sign, year, month, day, h, m, s, us)``
        * ``year``, ``month``, ``day``, ``hour``, ``minute``, ``second``
        * ``weekday`` — day of week as int (1 = Sunday … 7 = Saturday)
        * ``weekday_name`` — English name (``'sunday'``, ``'monday'``, …)

    Examples
    --------
    ISO string (proleptic Gregorian, BCE):

    >>> import montu
    >>> mtime = montu.Time('-2500-01-01 12:00:00.00')
    >>> print(mtime.jed)
    630823.0

    Julian Day Number:

    >>> mtime = montu.Time(807954.0, format='jd', scale='utc')

    Terrestrial-time seconds (J2000 epoch):

    >>> mtime = montu.Time(0, format='tt', scale='tt')
    >>> print(mtime.readable.datepro)
    '2000-01-01 11:58:56.0'

    Arithmetic — add ephemeris (TT) seconds with ``+`` / ``-``:

    >>> mtime3 = mtime + 365 * montu.DAY

    Arithmetic — shift by calendar units with :meth:`add` / :meth:`subs`:

    >>> mtime2 = mtime.add(years=1, days=10)
    >>> mtime4 = mtime.subs(days=30)

    Difference — ephemeris seconds with ``-``; calendar units with :meth:`diff`:

    >>> tt_seconds = later - earlier
    >>> cal_diff = later.diff(earlier)
    >>> cal_diff.years, cal_diff.to_days()
    """

    def __init__(self,
                 date=None,
                 format='iso',
                 scale='utc',
                 calendar='proleptic',
                 full=False,
                 zone=0):
        """Initialise a Time object.

        See class docstring for parameter and attribute descriptions.
        """
        # Representation is a dictionary with the representation
        self.readable = ReadableTime(self)

        # If date is None take now
        if date is None:
            date = datetime.now().isoformat().replace('T',' ')
            format = 'iso'
            scale = 'utc'
            calendar = 'proleptic'

        # If only a number is provided we assume is a terrestrial time
        if type(date) != str and format == 'iso':
            format = 'tt'
            scale = 'tt'

        """
        # Date is provided in Sothic (civil egyptian)
        if calendar == 'sothic':

            format = 'iso'
            calendar = 'mixed'
            pass
        """

        # Set calendar assumed
        self.calendar = calendar
        
        if format == 'iso':

            # Date was provided as a string
            if calendar=='proleptic':

                # Parse string
                self._parse_datestr(date)    

                # Calculate deltat: you only need year and month
                deltat = ROUND_SECONDS(
                    pymeeus_Epoch.tt2ut(
                    self.readable.comps[0]*self.readable.comps[1],
                    self.readable.comps[2])
                )

                # Proleptic Gregorian Julian day from parsed civil date
                jd = ROUND_JULIANDAYS(
                    _jed_from_datetime64(self.readable.obj_datetime64)
                )

                # According to scale, derive TT and UTC ephemeris seconds
                if scale == 'tt':
                    jtd = jd
                    tt = _et_from_jed(jtd)
                    et = tt - deltat
                    jed = ROUND_JULIANDAYS(jtd - deltat / DAY)
                else:
                    jed = jd
                    et = _et_from_jed(jed)
                    tt = et + deltat
                    jtd = ROUND_JULIANDAYS(jed + deltat / DAY)

                # Correct tt according to year interval
                year = self.readable.comps[0]*self.readable.comps[1]
                if (300<=year<=1582) and True:
                    jed_correction = JED_CORRECTION(year)
                    tt -= jed_correction
                    tt = ROUND_SECONDS(tt)

                # Initialize object according to tt
                self.update_time(tt,format='tt',scale='tt')

            elif calendar=='mixed':

                # Parse string
                self._parse_datestr(date)    

                # Calculated deltat
                deltat = ROUND_SECONDS(
                    pymeeus_Epoch.tt2ut(
                    self.readable.comps[0]*self.readable.comps[1],
                    self.readable.comps[2])
                )

                # Convert from components to pymeeus epoch which 
                # receives date in mixed calendar
                args = (self.readable.comps[0]*self.readable.comps[1],
                        self.readable.comps[2],
                        self.readable.comps[3]+\
                        self.readable.comps[4]/24+\
                        self.readable.comps[5]/(24*60)+\
                        self.readable.comps[6]/(24*60*60)+\
                        self.readable.comps[7]/(24*60*60*1e6))
                pymeeus_epoch = pymeeus_Epoch(*args)
                jd = pymeeus_epoch.jde()
                et = _et_from_jed(ROUND_JULIANDAYS(jd))
                
                # According to scale choose terrestrial time
                if scale == 'tt':
                    tt = et
                    jtd = ROUND_JULIANDAYS(jd)
                    et = tt - deltat
                    jed = ROUND_JULIANDAYS(jtd - deltat/DAY)
                else:
                    et = et
                    jed = ROUND_JULIANDAYS(jd)
                    tt = et + deltat
                    jtd = ROUND_JULIANDAYS(jed + deltat/DAY)

                # Initialize object according to tt
                self.update_time(tt,format='tt',scale='tt')

            elif calendar=='sothic':

                # Parse input date (optional trailing civil time)
                cdate = date.strip()
                djed = 0
                time_match = SOTHIC_TIME_RE.search(cdate)
                if time_match:
                    time_str = time_match.group(1)
                    cdate = cdate[:time_match.start()].strip()
                    hh, mm, ss = (float(c) for c in time_str.split(":"))
                    djed = (hh + mm/60 + ss/3600)/24

                # Compute horus day and julian day
                self.hed = Time._horus_date_to_days(cdate)
                self.jed = self.hed + JED_APOKATASTASIS + djed

                # Initialize object according to jed
                self.update_time(self.jed,format='jd',scale='utc')

            else:
                raise ValueError("Calendar '{calendar}' not recognzed. Use 'proleptic' or 'mixed'.")

            # Get readable
            self.get_readable()

        else:
            # Initialize object according to date and format and 
            self.update_time(date,format,scale)

            if full:
                self.get_readable()

        # Parse zone parameter
        zone_hours = 0.0
        if zone is not None and zone != 0:
            if isinstance(zone, montu.Observer) or hasattr(zone, 'lon'):
                zone_hours = zone.lon / 15.0
            elif isinstance(zone, str):
                if zone.upper().startswith('UTC'):
                    zone_val = zone[3:]
                else:
                    zone_val = zone
                
                # Check for colon format like UTC-5:30 or -5:30
                if ':' in zone_val:
                    parts = zone_val.split(':')
                    h = float(parts[0])
                    m = float(parts[1])
                    sign = -1.0 if '-' in parts[0] else 1.0
                    zone_hours = h + sign * m / 60.0
                else:
                    zone_hours = float(zone_val)
            else:
                zone_hours = float(zone)

        if zone_hours != 0:
            zonedt = zone_hours * montu.HOUR
            ft = self - zonedt
            ft.get_readable()
            self.__dict__.update(ft.__dict__)
            if hasattr(self, 'readable') and isinstance(self.readable, ReadableTime):
                object.__setattr__(self.readable, '_time_obj', self)

    def _parse_datestr(self,date):
        """Parse an ISO-like date string and populate ``self.readable``.

        Parameters
        ----------
        date : str
            Accepted formats::

                '2000-01-01'
                '2000-01-01 12:00:00'
                '2000-01-01 12:00:00.00'
                '-2500-01-01 12:00:00.00'   # proleptic BCE
                'bce2501-01-01 12:00:00.00' # calendar BCE
                '2501 bce 01-01 12:00:00.00'

        Notes
        -----
        Populates ``self.readable`` with ``datepro``, ``obj_datetime64``,
        ``comps``, ``year``, ``month``, ``day``, ``hour``, ``minute``,
        ``second``, ``usecond``, and ``datespice``.
        """

        # Default format
        style = 'astro' # Default style of input string 
        self.readable.datepro = date # Default date

        # Strip blank spaces
        date = date.strip()

        # Is time before current era
        self.bce = False
        if date[0] == '-':
            self.bce = True
        elif 'b' in date.lower():
            self.bce = True
            style = 'calendar'

        # Convert all formats to dateastro
        if self.bce and (style == 'calendar'):
            subs1 = lambda m:str(-(int(m.group(1))-1))
            subs2 = lambda m:str(-(int(m.group(1))-1))+'-'
            self.readable.datepro = re.sub(r'^bce[a-z\s]*(\d+)',subs1,self.readable.datepro.lower())
            self.readable.datepro = re.sub(r'(\d+)\s*b[\.]*c[\.]*e[\.]*\s*',subs2,self.readable.datepro.lower())

        # Create calendar and datetime object
        self.readable.obj_datetime64 = np.datetime64(self.readable.datepro)
        self.readable.comps = montu.Util.dt2cal(
            np.datetime64(self.readable.datepro.strip('-')),
            bce=self.bce)
        
        # Generate info
        self.readable.year = self.readable.comps[0]*self.readable.comps[1]
        self.readable.month = self.readable.comps[2]
        self.readable.day = self.readable.comps[3]
        self.readable.hour = self.readable.comps[4]
        self.readable.minute = self.readable.comps[5]
        self.readable.second = self.readable.comps[6]
        self.readable.usecond = self.readable.comps[7]
        
        # Adjust SPICE string according to epoch
        if self.bce:
            self.readable.datespice = f'{-self.readable.year+1:04d} B.C. {self.readable.month:02d}-{self.readable.day:02d} {self.readable.hour:02d}:{self.readable.minute:02d}:{self.readable.second:02d}.{self.readable.usecond:02d}'
        elif 0<self.readable.year<1000:
            self.readable.datespice = f'{self.readable.year:04d} A.D. {self.readable.month:02d}-{self.readable.day:02d} {self.readable.hour:02d}:{self.readable.minute:02d}:{self.readable.second:02d}.{self.readable.usecond:02d}'
        else:
            self.readable.datespice = self.readable.datepro
    
    def update_time(self,time=None,format='tt',scale='tt'):
        """Recompute all numeric time attributes from a single input value.

        This is the core internal method called by ``__init__``, arithmetic
        operators, and :meth:`get_readable`.

        Parameters
        ----------
        time : float, optional
            Numeric value in the scale given by *format*.  Defaults to
            ``self.tt`` (re-synchronise from current TT).
        format : {'tt', 'jd'}, optional
            Interpretation of *time*: ``'tt'`` for ephemeris seconds past
            J2000 (TT scale), ``'jd'`` for Julian Day number.
        scale : {'tt', 'utc'}, optional
            Whether *time* is in the TT or UTC scale.

        Notes
        -----
        Sets ``tt``, ``et``, ``jed``, ``jtd``, ``hed``, ``htd``,
        ``deltat``, ``bce``, ``isjulian``, ``obj_pyephem``, and
        ``obj_pyplanet``.
        """
        # Use if you set the attribute tt manually
        if time is None:
            time = self.tt
            format = 'tt'
            scale = 'tt'

        # Choose input format
        if format == 'jd':
            jd = time
        elif format == 'tt':
            et = time
            jd = _jed_from_et(et)
        else:
            raise AssertionError(f"Format '{format}' not recognized (valid 'iso', 'tt', 'jd')")

        # Initialize Epoch in pymeeus
        pymeeus_epoch = pymeeus_Epoch(jd)
        year,month,day = pymeeus_epoch.get_date()
        self.isjulian = pymeeus_Epoch.is_julian(year,month,day)
        self.deltat = ROUND_SECONDS(pymeeus_Epoch.tt2ut(year,month))
        self.bce = True if year<=0 else False
        
        # Get terrestrial time
        et = _et_from_jed(jd)
        if scale == 'tt':
            self.jtd = ROUND_JULIANDAYS(jd)
            self.tt = et
            self.jed = ROUND_JULIANDAYS(jd - self.deltat/DAY)
            self.et = self.tt - self.deltat
        else:
            self.jed = ROUND_JULIANDAYS(jd)
            self.et = et
            self.jtd = ROUND_JULIANDAYS(self.jed + self.deltat/DAY)
            self.tt = self.et + self.deltat

        # PyEphem Date: you need to provide jd with no deltat correction: is internal
        self.obj_pyephem = pyephem.Date(self.jed - PYEPHEM_JD_REF)

        # Create pyplanet epoch: you need to provide jd with no deltat correction: is internal
        self.obj_pyplanet = pyplanets_Epoch(self.jed)

        # Horus day (days since bce 2782-07-20 00:00:00)
        self.hed = self.jed - JED_APOKATASTASIS
        self.htd = self.jtd - JED_APOKATASTASIS

        self._set_weekday()
    
    def _set_weekday(self):
        """Set ``readable.weekday`` and ``readable.weekday_name`` from ``jed``."""
        index = int(self.jed + 1.5) % 7
        self.readable.weekday = index + 1
        self.readable.weekday_name = WEEKDAY_NAMES[index]

    @staticmethod
    def jed_to_weekday(jed, name=False):
        """Return the day of week for a Julian Day (UTC).

        Neither PyEphem nor PyMeeus expose weekday routines; this uses the
        standard integer-JD formula (same convention as the mixed calendar).

        Parameters
        ----------
        jed : float
            Julian Day in UTC.
        name : bool, optional
            If ``True``, return the English name (lowercase).  Otherwise
            return an integer with Sunday = 1 and Saturday = 7.

        Returns
        -------
        int or str

        Examples
        --------
        >>> import montu
        >>> montu.Time.jed_to_weekday(2451544.5)
        7
        >>> montu.Time.jed_to_weekday(2451545.5, name=True)
        'sunday'
        """
        index = int(jed + 1.5) % 7
        if name:
            return WEEKDAY_NAMES[index]
        return index + 1

    def get_readable(self):
        """Populate all human-readable string representations.

        Converts the current numeric time attributes to strings in all
        supported calendar systems and stores them in ``self.readable``.

        Returns
        -------
        self : Time
            Returns *self* so the call can be chained.

        Examples
        --------
        >>> import montu
        >>> mtime = montu.Time('-1000-03-21 06:00:00').get_readable()
        >>> print(mtime.readable.datepro)
        '-1000-03-21 06:00:00.0'
        >>> print(mtime.readable.datesot)
        '[hrw ...] ... ... ...'
        """

        # Update time
        self.update_time()

        # String for datemixed
        pyephem_str = f'{self.obj_pyephem}'.strip('-')
        parts = pyephem_str.split(' ')
        cals = [int(p) for p in parts[0].split('/')] +[int(p) for p in parts[1].split(':')] 
        if self.bce:
            # Adjust year if bce
            cals[0] -= 1
            cals[0] *= -1
        self.readable.datemix = f'{cals[0]}-{cals[1]:02d}-{cals[2]:02d} {cals[3]:02d}:{cals[4]:02d}:{cals[4]:02d}'

        # Human-readable civil date for the active calendar
        if self.calendar == 'proleptic':
            datestr = _datestr_from_datetime64(_datetime64_from_jed(self.jed))
        else:
            pymeeus_epoch = pymeeus_Epoch(self.jtd)
            year, month, day, hour, minute, second = pymeeus_epoch.get_full_date()
            datestr = (
                f'{year}-{int(month):02d}-{int(day):02d} '
                f'{int(hour):02d}:{int(minute):02d}:{second:04.1f}'
            )

        # Convert to sothic
        self.readable.datesot = Time._jed_to_sothic(self.jed)

        # Parse string
        self._parse_datestr(datestr)
        return self

    def __copy__(self):
        return Time(self.tt)

    def __add__(self, dtt):
        """Add ephemeris seconds (TT scale).

        Parameters
        ----------
        dtt : float
            Seconds to add in the TT scale (e.g. ``100 * montu.YEAR``).

        Returns
        -------
        Time
            New :class:`Time` advanced by *dtt* TT seconds.

        See Also
        --------
        add : shift by calendar years, months, days, etc.
        """
        if isinstance(dtt, Time):
            raise TypeError(
                "Use calendar units with Time.add() or ephemeris seconds "
                "with a numeric offset, not Time + Time."
            )
        new = copy.copy(self)
        new.tt += dtt
        new.update_time()
        return new

    def __radd__(self, dtt):
        return self.__add__(dtt)

    def __sub__(self, other):
        """Subtract ephemeris seconds or compute TT difference between instants.

        ``time - seconds`` moves *time* backward on the TT scale.
        ``time_a - time_b`` returns the TT difference in seconds (float).

        Parameters
        ----------
        other : float or Time
            Ephemeris seconds to subtract, or another :class:`Time`.

        Returns
        -------
        Time or float
            New :class:`Time` when *other* is numeric; TT seconds when both
            operands are :class:`Time`.

        See Also
        --------
        subs : shift backward by calendar units.
        diff : calendar-component difference.
        """
        if isinstance(other, Time):
            return self.tt - other.tt
        new = copy.copy(self)
        new.tt -= other
        new.update_time()
        return new

    def _calendar_kwargs(self, years=0, months=0, days=0, weeks=0,
                         hours=0, minutes=0, seconds=0):
        days += weeks * 7
        return dict(
            years=years, months=months, days=days,
            hours=hours, minutes=minutes, seconds=seconds,
        )

    def _shift_calendar(self, years=0, months=0, days=0, weeks=0,
                        hours=0, minutes=0, seconds=0):
        """Return a new Time shifted by calendar units."""
        days += weeks * 7
        calendar = getattr(self, 'calendar', 'proleptic')
        if calendar == 'sothic':
            return self._shift_sothic(
                years=years, months=months, days=days,
                hours=hours, minutes=minutes, seconds=seconds,
            )
        if calendar == 'mixed':
            return self._shift_mixed(
                years=years, months=months, days=days,
                hours=hours, minutes=minutes, seconds=seconds,
            )
        return self._shift_proleptic(
            years=years, months=months, days=days,
            hours=hours, minutes=minutes, seconds=seconds,
        )

    def _shift_proleptic(self, years=0, months=0, days=0, weeks=0,
                         hours=0, minutes=0, seconds=0):
        self._ensure_readable()
        era, y, m, d = self.readable.comps[:4]
        h, mi, s = self.readable.comps[4:7]
        astro_y = era * y + years
        m += months
        astro_y, m = _normalize_month_year(astro_y, m)
        y, m, d, h, mi, s = _apply_civil_offset(
            astro_y, m, d, h, mi, s,
            days=days, weeks=weeks,
            hours=hours, minutes=minutes, seconds=seconds,
            calendar='proleptic',
        )

        datepro = _mixed_datestr(y, m, d, h, mi, s)
        new = Time(datepro, calendar='proleptic')
        new.calendar = self.calendar
        return new

    def _shift_mixed(self, years=0, months=0, days=0, weeks=0,
                     hours=0, minutes=0, seconds=0):
        """Shift a mixed-calendar instant using civil months/years and UTC-day steps.

        Day-sized steps use Julian Day (``jed``) so the historical gap at the
        1582 Gregorian reform is respected (1582-10-04 + 1 day → 1582-10-15).
        """
        self._ensure_readable()
        jed = self.jed
        delta_days = days + weeks * 7 + hours / 24 + minutes / 1440 + seconds / 86400

        if years or months:
            y, m, d, h, mi, s = pymeeus_Epoch(jed).get_full_date()
            y, m, d, h, mi, s = _apply_civil_offset(
                y, m, d, h, mi, s,
                years=years, months=months, days=0, weeks=0,
                hours=0, minutes=0, seconds=0,
                calendar='mixed',
            )
            jed = pymeeus_Epoch(y, m, d, h, mi, s).jde()

        if delta_days:
            jed += delta_days

        return Time(jed, format='jd', scale='utc', calendar='mixed')

    def _shift_sothic(self, years=0, months=0, days=0,
                      hours=0, minutes=0, seconds=0):
        self._ensure_readable()
        total = self.hed + years * 365 + months * 30 + days
        total += (hours + minutes / 60 + seconds / 3600) / 24
        new = Time(total + JED_APOKATASTASIS, format='jd', scale='utc')
        new.calendar = 'sothic'
        return new

    def add(self, *, years=0, months=0, days=0, weeks=0,
            hours=0, minutes=0, seconds=0):
        """Advance *self* by calendar units.

        Parameters
        ----------
        years, months, days, weeks, hours, minutes, seconds : int or float, optional
            Calendar offsets in the active calendar system of *self*.

        Returns
        -------
        Time
            New :class:`Time` instance.

        Notes
        -----
        For uniform ephemeris (TT) offsets use ``+`` instead, e.g.
        ``mtime + 100 * montu.YEAR``.

        Examples
        --------
        >>> import montu
        >>> mtime = montu.Time('-1000-01-01 00:00:00')
        >>> mtime2 = mtime.add(years=1, days=10)
        """
        return self._shift_calendar(
            years=years, months=months, days=days, weeks=weeks,
            hours=hours, minutes=minutes, seconds=seconds,
        )

    def subs(self, *, years=0, months=0, days=0, weeks=0,
             hours=0, minutes=0, seconds=0):
        """Move *self* backward by calendar units.

        Parameters
        ----------
        years, months, days, weeks, hours, minutes, seconds : int or float, optional
            Calendar offsets in the active calendar system of *self*.

        Returns
        -------
        Time
            New :class:`Time` instance.

        Notes
        -----
        For uniform ephemeris (TT) offsets use ``-`` instead, e.g.
        ``mtime - montu.DAY``.

        Examples
        --------
        >>> import montu
        >>> mtime2 = mtime.subs(days=30)
        """
        return self._shift_calendar(
            years=-years, months=-months, days=-days, weeks=-weeks,
            hours=-hours, minutes=-minutes, seconds=-seconds,
        )

    def sub(self, dtt=None, **kwargs):
        """Deprecated alias — use :meth:`subs` or ``-`` for TT seconds."""
        if dtt is not None:
            if kwargs:
                raise TypeError("Pass either ephemeris seconds to sub() or calendar keywords, not both.")
            return self - dtt
        return self.subs(**kwargs)

    def _calendar_components(self):
        self._ensure_readable()
        calendar = getattr(self, 'calendar', 'proleptic')
        if calendar == 'sothic':
            hy, month, season, day = Time.parse_datesot(self.readable.datesot)
            total_days = Time._horus_days(
                hy,
                month,
                'MESUT' if season == 'mesut' else season.upper(),
                day,
            )
            day_fraction = self.jed - np.floor(self.jed)
            seconds = int(round(day_fraction * DAY))
            hours, rem = divmod(seconds, 3600)
            minutes, seconds = divmod(rem, 60)
            years, months, _, days = _sothic_civil_days_to_parts(total_days)
            return dict(
                years=years, months=months, days=days,
                hours=hours, minutes=minutes, seconds=seconds,
            )
        if calendar == 'mixed':
            y, m, d, h, mi, s = pymeeus_Epoch(self.jtd).get_full_date()
            day_int = int(d)
            day_fraction = d - day_int
            seconds = int(round((s + day_fraction * 86400)))
            hours, rem = divmod(seconds, 3600)
            minutes, seconds = divmod(rem, 60)
            return dict(
                years=int(y), months=int(m), days=day_int,
                hours=int(h) + hours, minutes=int(mi) + minutes, seconds=seconds,
            )
        era, y, m, d, h, mi, s, _ = self.readable.comps
        return dict(
            years=era * y, months=m, days=d,
            hours=h, minutes=mi, seconds=s,
        )

    @staticmethod
    def _component_diff(later_parts, earlier_parts, jed_days):
        years = later_parts['years'] - earlier_parts['years']
        months = later_parts['months'] - earlier_parts['months']
        days = later_parts['days'] - earlier_parts['days']
        hours = later_parts['hours'] - earlier_parts['hours']
        minutes = later_parts['minutes'] - earlier_parts['minutes']
        seconds = later_parts['seconds'] - earlier_parts['seconds']

        if seconds < 0:
            seconds += 60
            minutes -= 1
        if minutes < 0:
            minutes += 60
            hours -= 1
        if hours < 0:
            hours += 24
            days -= 1
        if days < 0:
            months -= 1
            ref_year = later_parts['years']
            ref_month = later_parts['months'] if months >= 0 else later_parts['months'] - 1
            if ref_month < 1:
                ref_month = 12
                ref_year -= 1
            days += _days_in_month(ref_year, ref_month)
        if months < 0:
            months += 12
            years -= 1

        return CalendarDelta(
            years=years, months=months, days=days,
            hours=hours, minutes=minutes, seconds=seconds,
            _jed_days=jed_days,
        )

    def diff(self, mtime):
        """Calendar difference ``self - mtime``.

        Parameters
        ----------
        mtime : Time
            Earlier or later instant to compare with *self*.

        Returns
        -------
        CalendarDelta
            Signed calendar components and Julian-day helpers.

        Notes
        -----
        For ephemeris (TT) seconds use ``self - mtime`` instead.

        Examples
        --------
        >>> import montu
        >>> t1 = montu.Time('2000-01-01')
        >>> t2 = montu.Time('2001-01-01')
        >>> t2.diff(t1).years
        1
        """
        if not isinstance(mtime, Time):
            raise TypeError("diff() expects another Time instance.")
        sign = 1
        later, earlier = self, mtime
        if self.tt < mtime.tt:
            later, earlier = mtime, self
            sign = -1
        jed_days = later.jed - earlier.jed
        delta = Time._component_diff(
            later._calendar_components(),
            earlier._calendar_components(),
            jed_days,
        )
        return -delta if sign < 0 else delta

    def _ensure_readable(self):
        """Populate human-readable fields when only partial data exist (e.g. weekday)."""
        if not hasattr(self.readable, "datepro"):
            self.get_readable()

    def details(self):
        """Return a string with the details of the Time object."""
        self._ensure_readable()

        str = f"""Montu Time Object:
-------------------------- 
Readable:
    Date in proleptic UTC (.readable.datepro): {self.readable.datepro}
    Date in mixed UTC (.readable.datemix): {self.readable.datemix}
    Date in SPICE format (.readable.datespice): {self.readable.datespice}
    Date in sothic format (.readable.datesot): {self.readable.datesot}
    Weekday (.readable.weekday): {self.readable.weekday} ({self.readable.weekday_name})
    Components (.readable.comps): {self.readable.comps}
Objects:
    Date in datetime64 format (.readable.obj_datetime64): {self.readable.obj_datetime64}
    Date in PyPlanet Epoch (.obj_pyplanet): {self.obj_pyplanet}
    Date in PyEphem Epoch (.obj_pyephem): {self.obj_pyephem}
General:
    Is bce (.bce): {self.bce}
    Is Julian (.isjulian): {self.isjulian}
Uniform scales:
    Terrestrial time:
        tt (.tt): {self.tt}
        jtd (.jtd): {self.jtd}
        htd (.htd): {self.htd}
    UTC time:
        et (.et): {self.et}
        jed (.jed): {self.jed}
        hed (.hed): {self.hed}
    Delta-t = TT - UTC (.deltat): {self.deltat}
"""
        print(str)        

    def __str__(self):

        self._ensure_readable()

        str = f"""{self.readable.datespice} / {self.readable.datesot}:
    Date in ISO format: {self.readable.datespice}
    Date in proleptic UTC: {self.readable.datepro}
    Date in mixed UTC: {self.readable.datemix}
    Weekday: {self.readable.weekday} ({self.readable.weekday_name})
    Date in sothic format: {self.readable.datesot}
    Terrestrial time: tt [seconds]: {self.tt}
    UTC time: jed [days]: {self.jed}
    Delta-t = TT - UTC [seconds]: {self.deltat}
"""
        return str
    
    def __repr__(self) -> str:
        self._ensure_readable()
        return (
            f"Time('{self.readable.datepro}'/'{self.readable.datemix}'/"
            f"'{self.readable.datesot}'/JED {self.jed}/JTD {self.jtd})"
        )

    def strftime(self,timefmt='%Y'):
        """Format the date as a string using ``strftime``-like codes.

        Parameters
        ----------
        timefmt : str, optional
            Format string with the following substitutions:

            * ``%Y`` - four-digit year (with sign for BCE)
            * ``%m`` - zero-padded month (01-12)
            * ``%d`` - zero-padded day (01-31)
            * ``%H`` - zero-padded hour (00-23)
            * ``%M`` - zero-padded minute (00-59)
            * ``%S`` - zero-padded second (00-59)

        Returns
        -------
        str
            Formatted date string.

        Examples
        --------
        >>> import montu
        >>> mtime = montu.Time('-1000-03-21 06:00:00').get_readable()
        >>> mtime.strftime('%Y-%m-%d')
        '-1000-03-21'
        """
        timefmt = timefmt.replace('%Y',f'{self.readable.year}')
        timefmt = timefmt.replace('%m',f'{self.readable.month:02d}')
        timefmt = timefmt.replace('%d',f'{self.readable.day:02d}')
        timefmt = timefmt.replace('%H',f'{self.readable.hour:02d}')
        timefmt = timefmt.replace('%M',f'{self.readable.minute:02d}')
        timefmt = timefmt.replace('%S',f'{self.readable.second:02d}')
        return timefmt
    
    @staticmethod
    def get_date(jed,format='comps'):
        """Return a date tuple or :class:`Time` object for a given Julian Day.

        Parameters
        ----------
        jed : float
            Julian Day (UTC scale).
        format : {'comps', 'mtime'}, optional
            ``'comps'`` (default) returns a raw PyEphem date tuple
            ``(year, month, day, hour, minute, second)``.
            ``'mtime'`` returns a :class:`Time` object.

        Returns
        -------
        tuple or Time
            Date components or :class:`Time` instance.

        Examples
        --------
        >>> import montu
        >>> montu.Time.get_date(2451545.0, format='comps')
        (2000, 1, 1, 12, 0, 0.0)
        >>> montu.Time.get_date(2451545.0, format='mtime')
        Time('2000-01-01 12:00:00'...)
        """
        comps = pyephem.Date(jed-montu.PYEPHEM_JD_REF).tuple()
        if format == 'comps':
            return comps
        elif format == 'mtime':
            datemix = f'{comps[0]}-{comps[1]:02d}-{comps[2]:02d} {comps[3]:02d}:{comps[4]:02d}:{comps[5]:06.3f}'
            mtime = montu.Time(datemix,calendar='mixed')
            return mtime
    
    @staticmethod
    def set_time_ticks(ax,format='tt',timefmt='%Y',**kwargs):
        """Label the x-axis ticks of a Matplotlib axes with formatted dates.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes whose x-tick positions will be relabelled.
        format : {'tt', 'jd'}, optional
            Interpret the tick values as TT seconds (``'tt'``) or Julian Days
            (``'jd'``).  Default ``'tt'``.
        timefmt : str, optional
            :meth:`strftime` format string for tick labels.  Default ``'%Y'``.
        **kwargs
            Extra keyword arguments forwarded to ``ax.set_xticklabels``.

        Examples
        --------
        >>> import matplotlib.pyplot as plt, montu
        >>> fig, ax = plt.subplots()
        >>> # (assume ax has been plotted with tt values on the x-axis)
        >>> montu.Time.set_time_ticks(ax, format='tt', timefmt='%Y-%m-%d')
        """
        tts = ax.get_xticks()
        xlabels = []
        for tt in tts:
            mtime = Time(tt,format=format).get_readable()
            xlabels += [f'{mtime.strftime(timefmt)}']
        ax.set_xticklabels(xlabels,**kwargs)

    @staticmethod
    def _horus_days(horus=0,month='I',season='Akhet',day=1):
        """Obtain the number of days since I-Akhet-1 of the
        year of the first ancient egyptian apokotástasis (2782 b.c.e.).

        Parameters:
            horus: int, default = 0:
                Horus year, i.e. civil year starting at 2782 bce.
            
            month: string, default = 'I':
                Month in the season. Values are I, II, III, IV

            season: string, default = 'Akhet':
                Season of the year. 
                Values are AKHET (A,AKH), PERET (P,PER), SHEMU (S,SHE).

            day: int, default = 1:
                Day. It must be between 1 and 30 (including).

        Return:
            Number of days since 2782 bce I-Akhet-1.
        """
        D = (int(day)-1) + (HORUS_SEASON[season]-1)*120 + (HORUS_MONTH[month]-1)*30 + int(horus)*365
        return D

    @staticmethod
    def _mixed_datemix_astronomical_year(datemix: str) -> int:
        """Return the astronomical year from a ``readable.datemix`` date part."""
        part = datemix.strip().split()[0]
        if part.startswith('-'):
            return -int(part.split('-')[1])
        return int(part.split('-')[0])

    @staticmethod
    def _mixed_year_tag_to_astronomical(year_tag: str) -> int:
        """Convert a sothic mixed-year tag to an astronomical calendar year."""
        tag = year_tag.strip()
        lower = tag.lower()
        if lower.startswith('bce'):
            year = int(re.sub(r'^bce\s*', '', lower, flags=re.I).strip())
            return -(year - 1)
        if lower.startswith('ce'):
            year = int(re.sub(r'^ce\s*', '', lower, flags=re.I).strip())
            return year
        return int(tag)

    @staticmethod
    def _horus_year_from_mixed_astronomical(astro_year: int) -> int:
        """Horus year whose civil I akhet 1 falls in the given mixed astronomical year."""
        lo = max(0, astro_year + 2700)
        hi = astro_year + 2850
        for horus in range(lo, hi):
            hd = Time._horus_days(horus, 'I', 'AKHET', 1)
            mtime = Time(hd + JED_APOKATASTASIS, format='jd', scale='utc', calendar='mixed')
            mtime.get_readable()
            if Time._mixed_datemix_astronomical_year(mtime.readable.datemix) == astro_year:
                return horus
        raise ValueError(
            f"No Horus year matches mixed calendar astronomical year {astro_year}."
        )

    @staticmethod
    def _parse_sothic_msd(sothic_date: str) -> tuple[int, str, str, int]:
        """Parse month, season, day and resolve Horus year from a sothic date string."""
        text = sothic_date.strip()
        match = SOTHIC_DATE_RE.match(text)
        if match:
            return (
                int(match.group(1)),
                match.group(2).upper(),
                match.group(3).upper(),
                int(match.group(4)),
            )

        match = SOTHIC_MIXED_DATE_RE.match(text)
        if match:
            astro_year = Time._mixed_year_tag_to_astronomical(match.group('year_tag'))
            horus = Time._horus_year_from_mixed_astronomical(astro_year)
            return (
                horus,
                match.group('month').upper(),
                match.group('season').upper(),
                int(match.group('day')),
            )

        raise ValueError(
            f"Cannot parse sothic date '{sothic_date}'. "
            "Expected '[hrw YEAR] MONTH SEASON DAY', "
            "'[bce YEAR] MONTH SEASON DAY', or '[-YEAR] MONTH SEASON DAY' "
            "(e.g. '[hrw 0] I akhet 1' or '[bce 2782] I akhet 1')."
        )

    @staticmethod
    def _horus_date_to_days(sothic_date):
        """Obtain the number of days since I-Akhet-1 of the
        year of the first ancient egyptian apokotástasis (2782 b.c.e.).

        Parameters:
            sothic_date: string:
                Horus-year format: ``[hrw HYEAR] MONTH SEASON DAY``

                or mixed-year format: ``[bce YEAR] MONTH SEASON DAY``,
                ``[ce YEAR] MONTH SEASON DAY``, ``[-YEAR] …``, ``[YEAR] …``
                where the bracketed year is the mixed Julian/Gregorian year
                of civil I akhet 1 for that Horus year.

                MONTH: I, II, III, IV
                SEASON: Akhet, Peret, Shemu (or Mesut for epagomenal days)
                DAY: 1–30 (1–5 for Mesut)

        Return:
            hd: int:
                Horus days. Number of days since 2782 bce I-Akhet-1.
        """
        horus, month, season, day = Time._parse_sothic_msd(sothic_date)

        if season == 'MESUT':
            if day < 1 or day > 5:
                raise ValueError(
                    f"Day '{day}' out of range for Mesut in '{sothic_date}'. "
                    "Epagomenal days must be between 1 and 5."
                )
            return horus * 365 + 360 + (day - 1)

        if month not in HORUS_MONTH.keys():
            raise ValueError(
                f"Month '{month}' not recognized in '{sothic_date}', "
                f"it must be among {tuple(HORUS_MONTH.keys())}"
            )
        if season not in HORUS_SEASON.keys():
            raise ValueError(
                f"Season '{season}' not recognized in '{sothic_date}', "
                f"it must be among {tuple(HORUS_SEASON.keys())} or Mesut"
            )
        if day < 1 or day > 30:
            raise ValueError(
                f"Day '{day}' out of range in '{sothic_date}'. "
                "It must be between 1 and 30."
            )

        return Time._horus_days(horus, month, season, day)

    @staticmethod
    def parse_datesot(datesot):
        """Parse ``readable.datesot`` into Horus year, month, season, and day.

        Returns season in lowercase (akhet, peret, shemu, mesut).
        """
        match = SOTHIC_DATE_RE.match(str(datesot).strip())
        if not match:
            raise ValueError(
                f"Cannot parse sothic date '{datesot}'. "
                "Expected format: '[hrw YEAR] MONTH SEASON DAY'."
            )
        hyear = int(match.group(1))
        month = match.group(2).upper()
        season = match.group(3).lower()
        day = int(match.group(4))
        return hyear, month, season, day

    @staticmethod
    def _horus_day_to_sothic(hd):
        """Convert from Horus day to sothic date

        Parameters:
            hd: int:
                Horus days, ie. number of days since
                2782 b.c.e. 07-20, which is the apokatastasis date
                according to the reference date of Censorino
                (see Lull, pag. 96)

        Return:
            cdate: string:
                Sothic date in the format ``[hrw HYEAR] MONTH SEASON DAY``
                where:
                    HYEAR: Horus year, 0 = 2782 bce.
                    MONTH: Number of month (I, II, III, IV)
                    SEASON: Name (letter of abreviation) of season:
                        Akhet (akh,a), Peret (per,p), Shemu (she,s)
                    DAY: Number of day
        """
        hy = int(hd/365)
        dy = hd%365
        if dy<360:
            s = int(dy/120) + 1
            ds = dy - 120*(s-1)
            m = int(ds/30) + 1
            d = int(ds - 30*(m-1)) + 1
            sothic = (
                f"[hrw {hy}] {MONTH_HORUS[str(m)]} "
                f"{SEASON_HORUS[str(s)].lower()} {d}"
            )
        else:
            # Epagomenos
            d = int(dy - 360 + 1)
            sothic = f"[hrw {hy}] I {SEASON_HORUS['0'].lower()} {d}"
        return sothic

    @staticmethod
    def _mixed_to_sothic(datemix):
        """Convert a mixed-calendar date string to a sothic date string.

        Parameters
        ----------
        datemix : str
            Date in the mixed calendar (``'CCYY-MM-DD HH:MM:SS'``).

        Returns
        -------
        str
            Egyptian civil date in the sothic format
            (e.g. ``'[hrw 300] II akhet 15'``).
        """
        # Create montu.Time object
        mtime = montu.Time(datemix,calendar='mixed')

        return Time._jed_to_sothic(mtime.jed)
    
    @staticmethod
    def _jed_to_sothic(jed):
        """Convert a Julian Day (UTC) to a sothic (Egyptian civil) date.

        Parameters
        ----------
        jed : float
            Julian Day in UTC.

        Returns
        -------
        str
            Egyptian civil date string (e.g. ``'[hrw 300] II akhet 15'``).

        Examples
        --------
        >>> montu.Time._jed_to_sothic(705497.5)
        '[hrw 0] I akhet 1'
        """
        # Compute yhe horus day
        hd = jed - JED_APOKATASTASIS

        # Compute the sothic date
        cdate = Time._horus_day_to_sothic(hd)
        return cdate

# (docstring moved into class body above)
