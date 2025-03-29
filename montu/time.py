
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
import spiceypy as spy
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

###############################################################
# Main class
###############################################################
class Time(object):

    def __init__(self,
                 date=None,
                 format='iso',
                 scale='utc',
                 calendar='proleptic',
                 full=False):
        
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

                # Convert to terrestrial time as if date was given in TT scale
                et = spy.utc2et(self.readable.datespice)
                deltat_leaps = spy.deltet(et,'ET')
                et -= deltat_leaps
                
                # Round et: we don't need precision below 1 second
                et = ROUND_SECONDS(et)
                
                # According to scale, add or remove deltat
                if scale == 'tt':
                    # et is terrestrial time
                    tt = et
                    et = et - deltat
                else:
                    # et is utc time
                    tt = et + deltat
                    et = et

                # Get Julian day
                jed = ROUND_JULIANDAYS(spy.unitim(et,'ET','JED'))
                jtd = ROUND_JULIANDAYS(spy.unitim(tt,'ET','JED'))

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
                et = ROUND_SECONDS(spy.unitim(jd,'JED','ET'))
                
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
        """Parse date string

        Parameters:
            date: string:
                Possible formats:
                    2000-01-01
                    2000-01-01 12:00:00
                    2000-01-01 12:00:00.00
                    -2500-01-01 12:00:00.00
                    bce2501-01-01 12:00:00.00
                    2501 bce 01-01 12:00:00.00 

        Update:
            Update the following attributes:
                bce: Is date before current era.
                datepro: Date in gregorian proleptic
                obj_datetime64: Object in datetime64
                components: Components of date
                year,month,day,hour,second,usecond: Components of date
                datespice: Date in SPICE format.
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
        """Update time object according to terrestrial time.
        
        This is the most important method.
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
            jd = ROUND_JULIANDAYS(spy.unitim(et,'ET','JED'))
        else:
            raise AssertionError(f"Format '{format}' not recognized (valid 'iso', 'tt', 'jd')")

        # Initialize Epoch in pymeeus
        pymeeus_epoch = pymeeus_Epoch(jd)
        year,month,day = pymeeus_epoch.get_date()
        self.isjulian = pymeeus_Epoch.is_julian(year,month,day)
        self.deltat = ROUND_SECONDS(pymeeus_Epoch.tt2ut(year,month))
        self.bce = True if year<=0 else False
        
        # Get terrestrial time
        et = ROUND_SECONDS(spy.unitim(jd,'JED','ET'))
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
    
    def get_readable(self):
        """Fill the readable attribute according to information.
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
        
        # Set string from terrestrial time
        self.readable.datespice = spy.et2utc(self.et+spy.deltet(self.et,'ET'),'C',4)
        datestr = self.readable.datespice

        # Converting from 
        sub_bc = lambda m:f'-{int(m.group(1))-1:04d}-{MONTH_ABREVS[m.group(2)]:02d}-'
        sub_ad = lambda m:f'{int(m.group(1)):04d}-{MONTH_ABREVS[m.group(2)]:02d}-'
        sub_nm = lambda m:f'-{MONTH_ABREVS[m.group(1)]:02d}-'
        if self.bce:
            datestr = re.sub(r'(\d+)\s*B.C.\s*(\w+)\s*',sub_bc,datestr)
        elif 'A.D.' in datestr:
            datestr = re.sub(r'(\d+)\s*A.D.\s*(\w+)\s*',sub_ad,datestr)
        else:
            datestr = re.sub(r'\s+(\w+)\s+',sub_nm,datestr)

        # Convert to caniucular
        self.readable.datecan = Time._jed_to_caniucular(self.jed)

        # Parse string
        self._parse_datestr(datestr)
        return self

    def __copy__(self):
        return Time(self.tt)

    def __add__(self,dtt):
        """Add terrestrial seconds without taking into accounts utc 
        corrections. This is the correct way to do it

        For example:
            mtime + 365*DAY -> will add 365 days worth-of-seconds to TT
        """
        new = copy.copy(self)
        new.tt += dtt
        new.update_time()
        return new
    
    def __sub__(self,dtt):
        new = copy.copy(self)
        new.tt -= dtt
        new.update_time()
        return new
    
    def add(self,dtt):
        """Add utc seconds. This is the way you should add
        seconds if want a given change in date.

        For example:
            mtime.add(365*DAY) -> will add 365 days to calendar
        """
        new = Time(self.jed + dtt/DAY,format='jd')
        return new
    
    def sub(self,dtt):
        """Substract utc seconds. This is the way you should substract
        seconds if want a given change in date.

        For example:
            mtime.add(365*DAY) -> will add 365 days to calendar
        """
        new = Time(self.jed - dtt/DAY,format='jd')
        return new
    
    def diff(self,mtime):
        """Compute difference between two dates in days
        """
        difference = self.jed - mtime.jed
        return difference

    def __str__(self):

        if len(self.readable.__dict__.keys()) == 0:
            self.get_readable()

        str = f"""Montu Time Object:
-------------------------- 
Readable:
    Date in proleptic UTC: {self.readable.datepro}
    Date in mixed UTC: {self.readable.datemix}
    Date in SPICE format: {self.readable.datespice}
    Date in caniucular format: {self.readable.datecan}
    Components: {self.readable.comps}
Objects:
    Date in datetime64 format: {self.readable.obj_datetime64}
    Date in PyPlanet Epoch: {self.obj_pyplanet}
    Date in PyEphem Epoch: {self.obj_pyephem}
General:
    Is bce: {self.bce}
    Is Julian: {self.isjulian}
Uniform scales:
    Terrestrial time:
        tt: {self.tt}
        jtd: {self.jtd}
        htd: {self.htd}
    UTC time:
        et: {self.et}
        jed: {self.jed}
        hed: {self.hed}
    Delta-t = TT - UTC = {self.deltat}
"""
        return str
    
    def __repr__(self) -> str:
        if len(self.readable.__dict__.keys()) == 0:
            str=f"Time(JED {self.jed}/JTD {self.jtd})"    
        else:
            str = f"Time('{self.readable.datepro}'/'{self.readable.datemix}'/'{self.readable.datecan}'/JED {self.jed}/JTD {self.jtd})"
        return str

    def strftime(self,timefmt='%Y'):
        timefmt = timefmt.replace('%Y',f'{self.readable.year}')
        timefmt = timefmt.replace('%m',f'{self.readable.month:02d}')
        timefmt = timefmt.replace('%d',f'{self.readable.day:02d}')
        timefmt = timefmt.replace('%H',f'{self.readable.hour:02d}')
        timefmt = timefmt.replace('%M',f'{self.readable.minute:02d}')
        timefmt = timefmt.replace('%S',f'{self.readable.second:02d}')
        return timefmt
    
    @staticmethod
    def get_date(jed,format='comps'):
        comps = pyephem.Date(jed-montu.PYEPHEM_JD_REF).tuple()
        if format == 'comps':
            return comps
        elif format == 'mtime':
            datemix = f'{comps[0]}-{comps[1]:02d}-{comps[2]:02d} {comps[3]:02d}:{comps[4]:02d}:{comps[5]:06.3f}'
            mtime = montu.Time(datemix,calendar='mixed')
            return mtime
    
    @staticmethod
    def set_time_ticks(ax,format='tt',timefmt='%Y',**kwargs):
        """Set xticks as Time objects
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
        comps[2] = comps[2].upper() # Season (akhet,peret,shemu)
        comps[3] = int(comps[3]) # Day (1..30)

        if comps[1] not in HORUS_MONTH.keys():
            raise ValueError(f"Month '{comps[1]}' not recognized in '{date}', it must be among {tuple(HORUS_MONTH.keys())}")
        if comps[2] not in HORUS_SEASON.keys():
            raise ValueError(f"Season '{comps[2]}' not recognized in '{date}', it must be among {tuple(HORUS_SEASON.keys())}")
        if (int(comps[3])>30) or (int(comps[3])<1):
            raise ValueError(f"Day '{comps[3]}' out of range. It must be between [1,30]")

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
        """Convert a mixed date to caniucular
        """
        # Create montu.Time object
        mtime = montu.Time(datemix,calendar='mixed')

        return Time._jed_to_caniucular(mtime.jed)
    
    @staticmethod
    def _jed_to_caniucular(jed):
        """Convert a mixed date to caniucular
        """
        # Compute yhe horus day
        hd = jed - JED_APOKATASTASIS

        # Compute the caniucular date
        cdate = Time._horus_day_to_caniucular(hd)
        return cdate

Time.__doc__ = """Create a time object
    
    This is one of the most important classes in the package, since
    it manages times and dates.

    Initialization parameters:
        date: string | float:

            Date provided. 
            
            When 'iso', possible formats are (all the same date):
                -1000-01-01 12:00:00.00
                bce1001-01-01 12:00:00.00
                1001 b.c.e. 01-01 12:00:00.00

            When 'sothic' (not implemented yet) the format is:
                20-II-shemu
            
            With seasons given by 5 leter names: shemu,akhet,peret

        format: string, default = 'iso':

            Format of the input date. 
            Possible values: 'iso', 'tt', 'jd', 'sothic'
        
        scale: string, default = 'utc':

            Scale of the time.
            Available: 'tt' (terrestrial time, uniform),  'utc' (coordinated, based on rotation).

        calendar: string, default = 'proleptic':

            Type of calendar of date. Proleptic gregorian correspond to the case when the Gregorian calendar
            is extended before the adoption date at 1582-10-15. When calendar = 'mixed' 
            the Julian calendar is used before the adoption date.

    Attributes:

        Time as strings:
        
            datepro: string:
                Date in gregorian prolectic, with format '[-]CCYY-MM-DD HH:MM:SS.fff'

            datemix: string:
                Date in gregorian mixed, with format '[-]CCYY-MM-DD HH:MM:SS.fff'

            datespice: string:
                Date in gregorian prolectic, with SPICE format.

        Time in uniform scales:
            deltat: float [seconds]
                Difference Dt = TT - UTC. 
            tt: float [seconds]
                Ephemeris time in tt. 
            et: float [seconds]
                Ephemeris time in utc. et = tt - deltat 
            jtd: float [days]
                Julian day in terrestrial time.
            jed: float [days]
                Julian day in UTC. jed = jtd - deltat/86400

        Time as special objects:

            obj_pyplanet: pyplanets.Epoch:
                Date in pyplanets format.

            obj_pyephem: pyephem.Date:
            
        datepro: string
            Date in gregorian prolectic but in astronomical format '[-]CCYY-MM-DD HH:MM:SS.fff'

        datemix: string
            Date in a mixed style (non-prolectic), with format '[bce]CCYY-MM-DD HH:MM:SS.fff', 
            meaning that when date is previous to 1582-10-15 the date is given in Julian 
            Calendar.

    Other attributes (related to frame of reference)

    Examples:

        Initialization using a string:

        
        Initialization using a float:

    """
