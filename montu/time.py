
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

# Egyptian civil calendar (caniucular)
HORUS_MONTH = dict(I=1,II=2,III=3,IV=4)
HORUS_SEASON = dict(AKHET=1,PERET=2,SHEMU=3,
                    A=1,P=2,S=3,
                    AKH=1,PER=2,SHE=3)
MONTH_HORUS = {str(v): k for k, v in HORUS_MONTH.items()}
SEASON_HORUS = {'0':'Mesut','1':'Akhet','2':'Peret','3':'Shemu'}

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

###############################################################
# Main class
###############################################################
class Time(object):
    """Represent an instant in time, supporting ancient astronomical calendars.

    This is the central class of MontuPython.  It converts dates between the
    proleptic Gregorian calendar, the mixed Julian/Gregorian calendar, the
    ancient Egyptian civil (caniucular) calendar, and several numerical
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

        When ``format='caniucular'``, the Egyptian civil date format is::

            '0-II-Akhet-1'                  # Year 0 (≈2782 BCE), I-Akhet-1

        When ``format='jd'``, *date* should be a Julian Day number (float).
        When ``format='tt'``, *date* should be ephemeris seconds past J2000.
    format : {'iso', 'tt', 'jd', 'caniucular'}, optional
        Format of *date*.  Default ``'iso'``.
    scale : {'utc', 'tt'}, optional
        Time scale of the input.  ``'utc'`` (default) treats *date* as a
        UTC epoch; ``'tt'`` treats it as Terrestrial Time (TT).
    calendar : {'proleptic', 'mixed', 'caniucular'}, optional
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
        * ``datecan`` — Egyptian civil (caniucular) date string
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

    Arithmetic — add calendar days:

    >>> mtime2 = mtime.add(365 * montu.DAY)

    Arithmetic — add TT seconds (no leap-second correction):

    >>> mtime3 = mtime + 365 * montu.DAY
    """

    def __init__(self,
                 date=None,
                 format='iso',
                 scale='utc',
                 calendar='proleptic',
                 full=False):
        """Initialise a Time object.

        See class docstring for parameter and attribute descriptions.
        """
        # Representation is a dictionary with the representation
        self.readable = montu.Dictobj()

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
        # Date is provided in Caniucular (civil egyptian)
        if calendar == 'caniucular':

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

            elif calendar=='caniucular':

                # Parse input date
                comps = date.split(' ')

                # Extract time
                if ':' in comps[-1]:
                    time = comps[-1]
                    hh,mm,ss = (float(c) for c in time.split(":"))
                    djed = (hh + mm/60 + ss/3600)/24
                else:
                    djed = 0

                # Check input format
                if '-' in comps[1].upper():
                    cdate = comps[0]+comps[1]
                else:
                    cdate = comps[0]

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
        >>> print(mtime.readable.datecan)
        'hrw ...-...-...-...'
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

        # Convert to caniucular
        self.readable.datecan = Time._jed_to_caniucular(self.jed)

        # Parse string
        self._parse_datestr(datestr)
        return self

    def __copy__(self):
        return Time(self.tt)

    def __add__(self,dtt):
        """Add a number of terrestrial seconds (TT scale, no leap corrections).

        Parameters
        ----------
        dtt : float
            Seconds to add in the TT scale.

        Returns
        -------
        Time
            New :class:`Time` object advanced by *dtt* seconds.

        Notes
        -----
        Use :meth:`add` instead if you want the result to shift the calendar
        date by exactly *dtt* seconds in UTC (accounting for ``deltat``).

        Examples
        --------
        >>> import montu
        >>> mtime = montu.Time('2000-01-01 00:00:00')
        >>> mtime2 = mtime + 365 * montu.DAY  # add 365 TT-days
        """
        new = copy.copy(self)
        new.tt += dtt
        new.update_time()
        return new
    
    def __sub__(self,dtt):
        """Subtract a number of terrestrial seconds (TT scale, no leap corrections).

        Parameters
        ----------
        dtt : float
            Seconds to subtract in the TT scale.

        Returns
        -------
        Time
            New :class:`Time` object moved back by *dtt* seconds.

        Examples
        --------
        >>> mtime2 = mtime - montu.DAY  # subtract one TT-day
        """
        new = copy.copy(self)
        new.tt -= dtt
        new.update_time()
        return new
    
    def add(self,dtt):
        """Add a number of UTC seconds, preserving calendar accuracy.

        Unlike the ``+`` operator, this method shifts the Julian Day in UTC
        and reconstructs the object, so the calendar date advances by the
        exact requested interval.

        Parameters
        ----------
        dtt : float
            Seconds to add in the UTC scale.

        Returns
        -------
        Time
            New :class:`Time` object advanced by *dtt* UTC seconds.

        Examples
        --------
        >>> import montu
        >>> mtime = montu.Time('-1000-01-01 00:00:00')
        >>> mtime2 = mtime.add(365 * montu.DAY)  # advance exactly one year
        """
        new = Time(self.jed + dtt/DAY,format='jd')
        return new
    
    def sub(self,dtt):
        """Subtract a number of UTC seconds, preserving calendar accuracy.

        Parameters
        ----------
        dtt : float
            Seconds to subtract in the UTC scale.

        Returns
        -------
        Time
            New :class:`Time` object moved back by *dtt* UTC seconds.

        Examples
        --------
        >>> import montu
        >>> mtime2 = mtime.sub(30 * montu.DAY)  # go back 30 calendar days
        """
        new = Time(self.jed - dtt/DAY,format='jd')
        return new
    
    def diff(self,mtime):
        """Return the Julian-day difference between two :class:`Time` objects.

        Parameters
        ----------
        mtime : Time
            The other :class:`Time` instance to subtract from *self*.

        Returns
        -------
        float
            ``self.jed - mtime.jed`` in fractional days.

        Examples
        --------
        >>> import montu
        >>> t1 = montu.Time('2000-01-01')
        >>> t2 = montu.Time('2000-03-01')
        >>> t2.diff(t1)
        60.0
        """
        difference = self.jed - mtime.jed
        return difference

    def _ensure_readable(self):
        """Populate human-readable fields when only partial data exist (e.g. weekday)."""
        if not hasattr(self.readable, "datepro"):
            self.get_readable()

    def __str__(self):

        self._ensure_readable()

        str = f"""Montu Time Object:
-------------------------- 
Readable:
    Date in proleptic UTC (.readable.datepro): {self.readable.datepro}
    Date in mixed UTC (.readable.datemix): {self.readable.datemix}
    Date in SPICE format (.readable.datespice): {self.readable.datespice}
    Date in caniucular format (.readable.datecan): {self.readable.datecan}
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
        return str
    
    def __repr__(self) -> str:
        self._ensure_readable()
        return (
            f"Time('{self.readable.datepro}'/'{self.readable.datemix}'/"
            f"'{self.readable.datecan}'/JED {self.jed}/JTD {self.jtd})"
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
    def _horus_date_to_days(caniucular_date):
        """Obtain the number of days since I-Akhet-1 of the
        year of the first ancient egyptian apokotástasis (2782 b.c.e.).

        Parameters:
            caniucular_date: string:
                Format: CCYY-MM-SS-DD

                where: 
                    CCYY: Horus year, 0 = 2782 bce.
                    MM: Number of month (I, II, III, IV)
                    SS: Name (letter of abreviation) of season:
                        Akhet (akh,a), Peret (per,p), Shemu (she,s)
                    DD: Number of day

        Return:
            hd: int:
                Horus days. Number of days since 2782 bce I-Akhet-1.
        """
        comps = caniucular_date.split('-')
        
        # Adjust ranges
        comps[0] = int(comps[0].strip('hHrw ')) # Horus year
        comps[1] = comps[1].upper() # Month (I, II, III, IV)
        comps[2] = comps[2].upper() # Season (Akhet, Peret, Shemu, Mesut)
        comps[3] = int(comps[3]) # Day (1..30, or 1..5 for Mesut)

        if comps[2] == 'MESUT':
            if comps[3] < 1 or comps[3] > 5:
                raise ValueError(
                    f"Day '{comps[3]}' out of range for Mesut in '{caniucular_date}'. "
                    "Epagomenal days must be between 1 and 5."
                )
            return int(comps[0]) * 365 + 360 + (comps[3] - 1)

        if comps[1] not in HORUS_MONTH.keys():
            raise ValueError(
                f"Month '{comps[1]}' not recognized in '{caniucular_date}', "
                f"it must be among {tuple(HORUS_MONTH.keys())}"
            )
        if comps[2] not in HORUS_SEASON.keys():
            raise ValueError(
                f"Season '{comps[2]}' not recognized in '{caniucular_date}', "
                f"it must be among {tuple(HORUS_SEASON.keys())} or Mesut"
            )
        if (int(comps[3])>30) or (int(comps[3])<1):
            raise ValueError(
                f"Day '{comps[3]}' out of range in '{caniucular_date}'. "
                "It must be between 1 and 30."
            )

        return Time._horus_days(*comps)

    @staticmethod
    def _horus_day_to_caniucular(hd):
        """Convert from Horus day to caniucular date

        Parameters:
            hd: int:
                Horus days, ie. number of days since
                2782 b.c.e. 07-20, which is the apokatastasis date
                according to the reference date of Censorino
                (see Lull, pag. 96)

        Return:
            cdate: string:
                Caniucular date in the format CCYY-MM-SS-DD where: 
                    CCYY: Horus year, 0 = 2782 bce.
                    MM: Number of month (I, II, III, IV)
                    SS: Name (letter of abreviation) of season:
                        Akhet (akh,a), Peret (per,p), Shemu (she,s)
                    DD: Number of day
        """
        prefix = 'hrw '
        hy = int(hd/365) 
        dy = hd%365
        if dy<360:
            s = int(dy/120) + 1
            ds = dy - 120*(s-1)
            m = int(ds/30) + 1
            d = int(ds - 30*(m-1)) + 1
            caniucular = f"{prefix}" + str(hy) + "-" + MONTH_HORUS[str(m)] + "-" + SEASON_HORUS[str(s)] + "-" + str(d)
        else:
            # Epagomenos
            d = int(dy - 360 + 1)
            caniucular = f"{prefix}" + str(hy) + "-" + "I" + "-" + SEASON_HORUS['0'] + "-" + str(d)
        return caniucular

    @staticmethod
    def _mixed_to_caniucular(datemix):
        """Convert a mixed-calendar date string to a caniucular date string.

        Parameters
        ----------
        datemix : str
            Date in the mixed calendar (``'CCYY-MM-DD HH:MM:SS'``).

        Returns
        -------
        str
            Egyptian civil date in the caniucular format
            (e.g. ``'hrw 300-II-Akhet-15'``).
        """
        # Create montu.Time object
        mtime = montu.Time(datemix,calendar='mixed')

        return Time._jed_to_caniucular(mtime.jed)
    
    @staticmethod
    def _jed_to_caniucular(jed):
        """Convert a Julian Day (UTC) to a caniucular (Egyptian civil) date.

        Parameters
        ----------
        jed : float
            Julian Day in UTC.

        Returns
        -------
        str
            Egyptian civil date string (e.g. ``'hrw 300-II-Akhet-15'``).

        Examples
        --------
        >>> montu.Time._jed_to_caniucular(705497.5)
        'hrw 0-I-Akhet-1'
        """
        # Compute yhe horus day
        hd = jed - JED_APOKATASTASIS

        # Compute the caniucular date
        cdate = Time._horus_day_to_caniucular(hd)
        return cdate

# (docstring moved into class body above)
