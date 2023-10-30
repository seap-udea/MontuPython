###############################################################
# Required packages
###############################################################
# Version 
from montu.version import *

# System packages
import os
import requests
import tqdm
import copy
import re

# Basic packages
import numpy as np
import spiceypy as spy
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from functools import lru_cache

# Special packages
from pyplanets.core.epoch import Epoch as pyplanets_Epoch
from pymeeus.Epoch import Epoch as pymeeus_Epoch
from pyplanets.core.coordinates import true_obliquity as pyplanets_true_obliquity
from pyplanets.core.coordinates import nutation_longitude as pyplanets_nutation_longitude
from pyplanets.core.coordinates import precession_equatorial
from pyplanets.core.angle import Angle as pyplanets_Angle

from pymeeus.Mercury import Mercury as pymeeus_Mercury
from pymeeus.Venus import Venus as pymeeus_Venus
from pymeeus.Earth import Earth as pymeeus_Earth
from pymeeus.Mars import Mars as pymeeus_Mars
from pymeeus.Jupiter import Jupiter as pymeeus_Jupiter
from pymeeus.Saturn import Saturn as pymeeus_Saturn
from pymeeus.Uranus import Uranus as pymeeus_Uranus
from pymeeus.Neptune import Neptune as pymeeus_Neptune

from pyplanets.planets.mercury import Mercury as pyplanets_Mercury
from pyplanets.planets.venus import Venus as pyplanets_Venus
from pyplanets.planets.earth import Earth as pyplanets_Earth
from pyplanets.planets.mars import Mars as pyplanets_Mars
from pyplanets.planets.jupiter import Jupiter as pyplanets_Jupiter
from pyplanets.planets.saturn import Saturn as pyplanets_Saturn
from pyplanets.planets.uranus import Uranus as pyplanets_Uranus
from pyplanets.planets.neptune import Neptune as pyplanets_Neptune

# Astropy
from astropy import constants as astropy_constants
from astropy import units as monu # Montu Units are astropy units
from astropy.time import Time as astropy_Time
from astroquery.jplhorizons import Horizons as astropy_Horizons
from astropy.coordinates import EarthLocation as astropy_EarthLocation
from astropy.coordinates import AltAz as astropy_AltAz
from astropy.coordinates import SkyCoord as astropy_SkyCoord
from astropy.coordinates import Angle as astropy_Angle

# Pyephem
import ephem as pyephem

# Avoid warnings
import warnings
warnings.filterwarnings("ignore")

###############################################################
# Global settings of the package
###############################################################
np.set_printoptions(precision=17)
pd.set_option("display.precision",17)

###############################################################
# Constants
###############################################################
# Styles: https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html#sphx-glr-gallery-style-sheets-style-sheets-reference-py
PLT_DEFAULT_STYLE = 'default' # others: ggplot, default, classic
SET_PLT_DEFAULT_STYLE = lambda:plt.style.use(PLT_DEFAULT_STYLE)

LEGACY = True
MAIN = False

# Numerical Constants
RAD = 180/np.pi
DEG = 1/RAD

# Phsyical
AU = 149597870.700 # km, value in Horizons
CSPEED = 299792.458 # km/s, value in Horizons

# Astronomical
MIN = 60 # s
HOUR = 60*MIN # s
DAY = 86400 # s
YEAR = 365.25*DAY # s
CENTURY = 100*YEAR # s
MILLENIUM = 10*YEAR # s

SIDEREAL_YEAR = 365.256363004*DAY # d, J2000
TROPICAL_YEAR = 365.242190402*DAY # d, J2000

# Rotation matrix between J2000 and ECLIPJ2000
M_J2000_ECLIPJ2000 = spy.pxform('J2000','ECLIPJ2000',0)

# Historical
PYEPHEM_JD_REF = 2415020.0
PYEPHEM_MJD_2000 = 36525.0
JED_2000 = 2451545.0

# Datum of the calendar in Julian year
"""Note from wikipedia:
'The Julian day number (JDN) is the integer assigned to a whole solar day 
in the Julian day count starting from noon Universal Time, with Julian day 
number 0 assigned to the day starting at noon on Monday, January 1, 4713 BC, 
proleptic Julian calendar (November 24, 4714 BC, in the proleptic Gregorian calendar)'
"""
JULIAN_DATUM_PROLEPTIC = 4713
JULIAN_DATUM_JULIAN = 4712

# Required kernels
"""This dictionary describe the kernels the package require to compute planetary positions

If the dictionary is blank it means that the kernel is in the data directory.
"""
KERNELS = {
    'latest_leapseconds.tls':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls',
    'de441_part-1.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp',
    'de441_part-2.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp',
    'frame.tk':'',
    'pck00011.tpc':'',
    'earth_assoc_itrf93.tf':''
}
KERNELS_LOADED = dict()

# Planets
PLANETARY_IDS = dict(
    SUN = 10,
    MERCURY = 1,
    VERNUS = 2,
    EARTH = 399,
    MOON = 301,
    MARS = 4,
    JUPITER = 5,
    SATURN = 6,
    URANUS = 7,
    NEPTUNE = 8,
)
PLANETARY_NAMES = {str(v): k for k, v in PLANETARY_IDS.items()}

###############################################################
# Montu Python Util Class
###############################################################
class Montu(object):

    def vprint(verbose,*args):
        """Print messages in verbose mode
        """
        if verbose:
            print(*args)

    def _data_path(filename,check=False):
        """Get the full path of the `datafile` which is one of the datafiles provided with the package.
        
        Parameters:
            filename: Name of the data file, string.
            
        Return:
            Full path to package datafile in the python environment.
            
        """
        file_path = os.path.join(os.path.dirname(__file__),'data',filename)
        if check and (not os.path.isfile(file_path)):
            raise ValueError(f"File '{filename}' does not exist in data directory")
        return file_path

    def _wget(url, filename, verbose=False):
        """Get a file from a url and store it as filename

        Source: ChatGPT
        """
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        Montu.vprint(verbose,f"Downloading {filename} from {url} [size = {total_size}]")

        # Initialize the progress bar
        progress_bar = tqdm.tqdm(total=total_size, unit='B', unit_scale=True)

        with open(filename, 'wb') as file:
            for data in response.iter_content(chunk_size=1024):
                file.write(data)
                progress_bar.update(len(data))

        progress_bar.close()

    def dt2cal(dt,bce=False):
        """Convert array of datetime64 to a calendar array of year, month, day, hour,
        minute, seconds, microsecond with these quantites indexed on the last axis.

        Parameters
            dt : datetime64 array (...)
                numpy.ndarray of datetimes of arbitrary shape

        Returns
            cal : uint32 array (..., 7)
            
                calendar array with last axis representing year, month, day, hour,
                minute, second, microsecond

        Adapted from: https://stackoverflow.com/a/56260054
        """
        out = np.empty(dt.shape + (7,), dtype="u4")
        Y, M, D, h, m, s = [dt.astype(f"M8[{x}]") for x in "YMDhms"]
        out[..., 0] = Y + 1970 # Gregorian Year
        
        out[..., 1] = (M - Y) + 1 # month
        out[..., 2] = (D - M) + 1 # dat
        out[..., 3] = (dt - D).astype("m8[h]") # hour
        out[..., 4] = (dt - h).astype("m8[m]") # minute
        out[..., 5] = (dt - m).astype("m8[s]") # second
        out[..., 6] = (dt - s).astype("m8[us]") # microsecond

        #out = np.array([float(o) for o in out])
        out = [int(o) for o in out]
        out = [-1]+out if bce else [1]+out
        return out

    def load_kernels(kernels=KERNELS,dir='/tmp/',verbose=True):
        
        # Check if dir exists
        if not os.path.exists(dir):
            os.system(f"mkdir -p {dir}")    

        # Load kernel
        for kernel,item in kernels.items():

            #
            if kernel in KERNELS_LOADED.keys():
                Montu.vprint(verbose,f"Kernel {kernel} already loaded, skipping")
                continue

            # Local kernel
            if len(item) == 0:
                if verbose:print(f"Loading local kernel {kernel}")
                kernel_file = Montu._data_path(kernel)
                if os.path.isfile(kernel_file):
                    spy.furnsh(kernel_file)
                else:
                    raise AssertionError(f"Kernel file '{kernel}' not found in data directory")
                KERNELS_LOADED[kernel] = True
                continue

            # Remote kernel
            kernel_path = dir+"/"+kernel
            
            if not os.path.exists(kernel_path):
                # Download kernel if it is not yet downloaded
                if verbose:print(f"Downloading '{kernel}'...")
                Montu._wget(item,kernel_path)
            
            # Once downloaded furnish kernel
            if verbose:print(f"Loading kernel {kernel}")
            spy.furnsh(kernel_path)
            KERNELS_LOADED[kernel] = True

    def print_df(df):
        """Print DataFrame.
        
        Parameters:
            df: Pandas DataFrame:
                DataFrame to print.
        """
        from IPython.display import display,HTML
        display(HTML(df.to_html()))

    def _linear_map(mapped,observed):
        a = (observed[1]-observed[0])/(mapped[1]-mapped[0])
        b = observed[0] - a*mapped[0]
        map = lambda x:a*x+b
        return map
    
    def dec2hex(dec,string=True):
        sgn = np.sign(dec)
        dec = abs(dec)
        h = int(dec)
        mf = 60*(dec - int(dec))
        m = int(mf)
        s = 60*(mf - m)
        if string:
            ret = f"{int(sgn*h):02d}:{int(m):02d}:{s:.3f}"
        else:
            ret = sgn*h,m,s
        return ret
    
    def string_difference(string1, string2):
        """Calculate the difference between two strings
        """
        A = set(string1.split()) 
        B = set(string2.split()) 
        str_diff = A.symmetric_difference(B)
        isEmpty = (len(str_diff) == 0)
        return str_diff
    
    def haversine_distance(lat1, lon1, lat2, lon2):
        """Compute angular distance between two points
        """
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        return c
    
    def montu_mark(ax):
        """Add a water mark to a 2d or 3d plot.
        
        Parameters:
        
            ax: Class axes: 
                Axe where the watermark will be placed.
        """
        #Get the height of axe
        axh=ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted()).height
        fig_factor=axh/4
        
        #Options of the water mark
        args=dict(
            rotation=270,ha='left',va='top',
            transform=ax.transAxes,color='pink',fontsize=6*fig_factor,zorder=100
        )
        
        #Text of the water mark
        mark=f"MontuPython {version}"
        
        #Choose the according to the fact it is a 2d or 3d plot
        try:
            ax.add_collection3d
            plt_text=ax.text2D
        except:
            plt_text=ax.text
            
        text=plt_text(1,1,mark,**args);
        return text

# Aliases
D2H = Montu.dec2hex

###############################################################
# Configuration class
###############################################################
class MontuConfiguration(object):
    Verbosity = False

###############################################################
# MonTime Class
###############################################################
class MonTime(object):
    """Create a time-frame object
    
    This is one of the most important classes in the package, since
    it manages times and dates.

    Initialization parameters:
        date: string | float
            Date. When iso, possible formats are (all the same date):
                -1000-01-01 12:00:00.00
                bce1001-01-01 12:00:00.00
                1001 b.c.e. 01-01 12:00:00.00

        format: string, default = 'iso'
            Format of the input date. Possible values: 'iso', 'tt', 'jd'
        
        scale: string, default = 'utc'
            Scale of the time: 'tt' (terrestrial time, uniform), 
            'utc' (coordinated, based on rotation).

        proleptic: boolean, default = True
            Proleptic gregorian correspond to the case when the Gregorian calendar
            is extended before the adoption date at 1582-10-15. 

            When proleptic = False, the Julian calendar is used before the adoption
            date.

    Attributes:

        Time as strings:
        
            datestr: string
                Date in gregorian prolectic, with format '[-]CCYY-MM-DD HH:MM:SS.fff'
                
        Time in uniform scales:
            tt: float [seconds]
                Ephemeris time in tt. 
            et: float [seconds]
                Ephemeris time in utc. 
            deltat: float [seconds]
                Difference Dt = TT - UTC. 
            jtd: float [days]
                Julian day in terrestrial time.
            jed: float [days]
                Julian day in UTC.

        Time as special objects:

            dateastro: string
                Date in astropy format '[-]CCYY-MM-DD HH:MM:SS.fff'
            
        datemix: string
            Date in a mixed style (non-prolectic), with format '[bce]CCYY-MM-DD HH:MM:SS.fff', 
            meaning that when date is previous to 1582-10-15 the date is given in Julian 
            Calendar.
        dateastro: string
            Date in gregorian prolectic but in astronomical format '[-]CCYY-MM-DD HH:MM:SS.fff'

    Other attributes (related to frame of reference)

        epsilon: float [deg]
            True obliquity for the date.
            
        delta_psi: float [deg]
            Nutation longitude for the date

        M_J2000_Epoch: float
            Matrix transforming from equinox J2000 to equinos at Epoch.
            
    Examples:

        Initialization using a string:

        
        Initialization using a float:

    """    
    def __init__(self,
                 date,
                 format='iso',
                 scale='utc',
                 calendar='proleptic',
                 verbose=False):
        
        if type(date) != str and format == 'iso':
            # If only a number is provided we assume is a terrestrial time
            format = 'tt'
            scale = 'tt'

        # Set calendar
        self.calendar = calendar
        
        if format == 'iso':

            # Date was provided as a string
            if calendar=='proleptic':
                # Parse string
                self._parse_datestr(date)

                # Calculated deltat
                deltat = pymeeus_Epoch.tt2ut(self.components[0]*self.components[1],
                                                    self.components[2])

                # Convert to TT
                et = spy.str2et(self.datespice)
                et -= round(spy.deltet(et,'ET'),3)
                et = round(et,3) # Precision in time is milliseconds
                
                if scale == 'tt':
                    tt = et
                    et = et - deltat
                else:
                    tt = et + deltat
                    et = et
                
                # Get Julian day
                jed = spy.unitim(et,'ET','JED')
                jtd = spy.unitim(tt,'ET','JED')

            elif calendar=='mixed':
                # Parse string
                self._parse_datestr(date)

                # Calculated deltat
                deltat = pymeeus_Epoch.tt2ut(self.components[0]*self.components[1],
                                                    self.components[2])

                # Create pyplanet object since input format is non-prolectic
                args = (self.components[0]*self.components[1],
                        self.components[2],
                        self.components[3]+\
                        self.components[4]/24+\
                        self.components[5]/(24*60)+\
                        self.components[6]/(24*60*60)+\
                        self.components[7]/(24*60*60*1e6))
                pyplanet = pyplanets_Epoch(*args)
                
                # According to scale
                jd = pyplanet.jde()
                et = spy.unitim(jd,'JED','ET')
                
                # Choose according to scale
                if scale == 'tt':
                    tt = et
                    jtd = jd
                    et = tt - deltat
                    jed = jtd - deltat/DAY
                else:
                    et = et
                    jed = jd
                    tt = et + deltat
                    jtd = jed + deltat/DAY
            else:
                raise ValueError("Calendar '{calendar}' not recognzed. Use 'proleptic', 'mixed'.")

            # Initialize object according to tt
            self.update_time(tt,format='tt',scale='tt')

        else:
            # Initialize object according to date and format and 
            self.update_time(date,format,scale)

    def _parse_datestr(self,date):
        """Parse string
        """
        # Default format
        style = 'astro' # Default style of input string 
        self.datestr = date # Default date

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
            self.datestr = re.sub('bce[a-z\s]*(\d+)',subs1,self.datestr.lower())
            self.datestr = re.sub('(\d+)\s*b[\.]*c[\.]*e[\.]*\s*',subs2,self.datestr.lower())

        # Create calendar and datetime object
        self.obj_datetime64 = np.datetime64(self.datestr)
        self.components = Montu.dt2cal(np.datetime64(self.datestr.strip('-')),
                                       bce=self.bce)
        # Pay attention to leap seconds
        try:
            self.obj_datetime = datetime(*self.components[1:])
        except ValueError:
            components = copy.deepcopy(self.components)
            components[3] = components[3] - 1
            self.obj_datetime = datetime(*components[1:])
        
        # Generate info
        self.year = self.components[0]*self.components[1]
        self.month = self.components[2]
        self.day = self.components[3]
        self.hour = self.components[4]
        self.minute = self.components[5]
        self.second = self.components[6]
        self.usecond = self.components[7]
        
        # Datetime in Julian years
        julian_year = (JULIAN_DATUM_PROLEPTIC - self.components[1]) + 1
        if julian_year>0:
            self.obj_datejulian = datetime(*((julian_year,)+tuple(self.components[2:])))
        else:
            julian_year -= 1 # To avoid year 0
            self.obj_datejulian = datetime(*((-julian_year,)+tuple(self.components[2:])))

        # Date in julian years
        self.datejulian = self.obj_datejulian.strftime('%Y-%m-%d %H:%M:%S.%f')
        
        # Convert to SPICE format
        if self.bce:
            datelist = [int(f) for f in self.obj_datetime.strftime('%Y,%m,%d,%H,%M,%S,%f').split(',')]
            datelist[0] += 1
            dateprev = datetime(*tuple(datelist))
        else:
            dateprev = self.obj_datetime

        # Adjust SPICE string according to epoch
        if self.bce:
            self.datespice = dateprev.strftime('%Y B.C. %m-%d %H:%M:%S.%f')
        elif 0<self.components[1]<1000:
            self.datespice = dateprev.strftime('%Y A.D. %m-%d %H:%M:%S.%f')
        else:
            self.datespice = self.datestr

    def update_time(self,time=None,format='tt',scale='tt'):
        """This is the most important method.
        """

        # Use if you set the attribute tt manually
        if time is None:
            time = self.tt
            format = 'tt'
            scale = 'tt'

        self.scale = scale
        # Initialize dates using terrestrial time
        if format == 'jd':
            jd = time
        elif format == 'tt':
            et = time
            jd = spy.unitim(et,'ET','JED')
        else:
            raise AssertionError(f"Format '{format}' not recognized (valid 'iso', 'tt', 'jd')")

        # Initialize Epoch
        pyplanet = pyplanets_Epoch(jd)
        year,month,day = pyplanet.get_date()
        self.deltat = pymeeus_Epoch.tt2ut(year,month)
        self.bce = True if year<=0 else False
        
        # Terrestrial time
        et = spy.unitim(jd,'JED','ET')
        if scale == 'tt':
            self.jtd = jd
            self.tt = et
            self.jed = jd - self.deltat/DAY
            self.et = self.tt - self.deltat
        else:
            self.jed = jd
            self.et = et
            self.jtd = self.jed + self.deltat/DAY
            self.tt = self.et + self.deltat

        # Create pyplanet epoch: you need to provide jd with no deltat correction: is internal
        self.obj_pyplanet = pyplanets_Epoch(self.jed)

        # Create astrotime 
        self.obj_astrotime = astropy_Time(self.jtd,format='jd',scale='tdb')
        
        # PyEphem Date
        self.obj_pyephem = pyephem.Date(self.jed - PYEPHEM_JD_REF)
        
        # Generate object and string for datemixed 
        pyephem_str = f'{self.obj_pyephem}'.strip('-')
        parts = pyephem_str.split(' ')
        cals = [int(p) for p in parts[0].split('/')] +[int(p) for p in parts[1].split(':')]

        # Mixed normal
        if self.bce:
            cals[0] = cals[0]-1
        # Be careful with the leap months
        try:
            self.obj_datetimemix = datetime(*cals)
        except ValueError:
            cals[1] = cals[1]-1
            self.obj_datetimemix = datetime(*cals)
        if self.bce:
            self.datemixed = self.obj_datetimemix.strftime('-%Y-%m-%d %H:%M:%S')
        else:
            self.datemixed = self.obj_datetimemix.strftime('%Y-%m-%d %H:%M:%S')
        
        # Mixed in Julian years
        julian_year = (JULIAN_DATUM_JULIAN - cals[0] - 1) + 1
        if julian_year>0:
            self.obj_datejulianmix = datetime(*((julian_year,)+tuple(cals[1:])))
        else:
            julian_year -= 1 # To avoid year 0
            self.obj_datejulianmix = datetime(*((-julian_year,)+tuple(cals[1:])))
        self.datejulianmixed = self.obj_datejulianmix.strftime('%Y-%m-%d %H:%M:%S')
        
        # Set string
        c = spy.et2utc(self.et+spy.deltet(self.et,'ET'),'C',4)
        sub = lambda m:f'{int(m.group(1)):04d} '
        if self.bce:
            c = re.sub('(\d+)\s*B.C.\s*',sub,c)
            # Be careful with leap year
            try:
                d = datetime.strptime(c,'%Y %b %d %H:%M:%S.%f')
            except ValueError:
                d = self.obj_datetimemix
            datestr = f'{-(d.year-1)}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}.{d.microsecond}'
        else:
            if 'A.D.' in c:
                c = re.sub('(\d+)\s*A.D.\s*',sub,c)
                d = datetime.strptime(c,'%Y %b %d %H:%M:%S.%f')
            else:
                d = datetime.strptime(c,'%Y %b %d %H:%M:%S.%f')
                
            datestr = f'{d.year}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}.{d.microsecond}'
        
        # Parse string
        self._parse_datestr(datestr)

        # True obliquity and nutation longitude 
        self.epsilon = float(pyplanets_true_obliquity(self.obj_pyplanet))
        self.delta_psi = float(pyplanets_nutation_longitude(self.obj_pyplanet))
        self.M_equatorial_ecliptic = spy.rotate(self.epsilon*DEG,1)

        # Greenwich True Sidereal Time (GTST)
        self.gtst = 24*self.obj_pyplanet.apparent_sidereal_time(self.epsilon,
                                                                self.delta_psi)

        # Update matrices
        self.M_J2000_Epoch = spy.pxform('J2000','EARTHTRUEEPOCH',self.tt)
        self.M_Epoch_J2000 = spy.invert(self.M_J2000_Epoch)
        self.M_EJ2000_Epoch = spy.pxform('ECLIPJ2000','EARTHTRUEEPOCH',self.tt)
        self.M_Epoch_EJ2000 = spy.invert(self.M_J2000_Epoch)

    def _get_signature(self):
        signature = ''
        keys = sorted(self.__dict__.keys())[::-1]
        for key in keys:
            item = self.__dict__[key]
            signature += key+':'+str(item)+'__'
        return signature

    def _get_hash(self):
        return str(hash(self._get_signature()))
    
    def is_different(self,mtime):
        return Montu.string_difference(self._get_signature(),mtime._get_signature())

    def __add__(self,dtt):
        new = copy.deepcopy(self)
        new.tt += dtt
        new.update_time()
        return new
    
    def __sub__(self,dtt):
        new = copy.deepcopy(self)
        new.tt -= dtt
        new.update_time()
        return new

    """
    def __sub__(self,mtime):
        return self.jed-mtime.jed
    """
    
    def __str__(self):
        str = f"""Montu Time Object:
--------------------------
General:
    Calendar: {self.calendar}
    Is bce: {self.bce}
    Components UTC: {self.components}
Uniform scales:
    Terrestrial time:
        tt: {self.tt}
        jtd: {self.jtd}
    UTC time:
        et: {self.et}
        jed: {self.jed}
    Delta-t = TT - UTC = {self.deltat}
Strings:
    Date in SPICE format: {self.datespice}
    Date in proleptic calendar: {self.datestr}
    Date in proleptic calendar (jul.year): {self.datejulian}
    Date in mixed calendar: {self.datemixed}
    Date in mixed calendar (jul.year): {self.datejulianmixed}
Objects:
    Date in datetime64 format: {self.obj_datetime64}
    Date in datetime format proleptic: {self.obj_datetime}
    Date in datetime format proleptic (julian year): {self.obj_datejulian}
    Date in datetime format mixed: {self.obj_datetimemix}
    Date in datetime format mixed (julian year): {self.obj_datejulianmix}
    Date in PyPlanet Epoch: {self.obj_pyplanet}
    Date in PyEphem Epoch: {self.obj_pyephem}
    Date in AstroPy Time: {self.obj_astrotime}
Astronomical properties at Epoch:
    True obliquity of ecliptic: {D2H(self.epsilon,1)}
    True nutation longitude: {D2H(self.delta_psi,1)}
    Greenwhich Meridian Sidereal Time: {D2H(self.gtst,1)}
Hash: {self._get_hash()}
"""
        return str
    
    def get_equinoxes_solstices_year(self):
        firstday_jed = pymeeus_Epoch(self.year,1,1).jde()
        firstdate = pyephem.Date(firstday_jed - PYEPHEM_JD_REF)
        vernal_jed = pyephem.next_vernal_equinox(firstdate) + PYEPHEM_JD_REF
        summer_jed = pyephem.next_summer_solstice(firstdate) + PYEPHEM_JD_REF
        auttumnal_jed = pyephem.next_autumnal_equinox(firstdate) + PYEPHEM_JD_REF
        winter_jed = pyephem.next_winter_solstice(firstdate) + PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed
    
    @staticmethod
    def next_equinoxes_solstices(mtime):
        date = pyephem.Date(mtime.jed - PYEPHEM_JD_REF)
        vernal_jed = pyephem.next_vernal_equinox(date) + PYEPHEM_JD_REF
        summer_jed = pyephem.next_summer_solstice(date) + PYEPHEM_JD_REF
        auttumnal_jed = pyephem.next_autumnal_equinox(date) + PYEPHEM_JD_REF
        winter_jed = pyephem.next_winter_solstice(date) + PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed
    
    @staticmethod
    def previous_equinoxes_solstices(mtime):
        date = pyephem.Date(mtime.jed - PYEPHEM_JD_REF)
        vernal_jed = pyephem.previous_vernal_equinox(date) + PYEPHEM_JD_REF
        summer_jed = pyephem.previous_summer_solstice(date) + PYEPHEM_JD_REF
        auttumnal_jed = pyephem.previous_autumnal_equinox(date) + PYEPHEM_JD_REF
        winter_jed = pyephem.previous_winter_solstice(date) + PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed
    
###############################################################
# Stars Class
###############################################################
class Stars(object):

    def __init__(self,data=None,filename=None):

        if data is not None:
            # Load data for stars from a dataframe already loaded
            self.data = copy.deepcopy(data)
            
        elif filename:
            # Load data from a file
            self.data = pd.read_csv(filename)

        else:
            # Load data from the database provided with package
            self.data = pd.read_csv(Montu._data_path('stars.csv',check=True))
        
        self.number = len(self.data)

    def get_stars(self,**args):
        """Filter stars by criteria

        Examples:
            # Get a single stars
            aldebaran = allstars.get_stars(ProperName='Aldebaran')

            # All visible stars in the sky
            visible = allstars.get_stars(Mag=[-1,6.5])

            # All visible stars with declination less than 1 deg
            equator = allstars.get_stars(Mag=[-1,6.5],Dec=[-1,1])
        """

        # If no args get all stars in data base
        if len(args)==0:
            return self.data
        
        # If args provided it will try to filter database according to conditions
        cond = np.array([True]*len(self.data))
        for key,item in args.items():
            if isinstance(item,list):
                min = float(item[0])
                max = float(item[1])
                cond = (self.data[key]>=min)&(self.data[key]<=max)&(cond)
            else:   
                cond = (self.data[key]==item)&(cond)
    
        return Stars(self.data[cond])
    
    def get_stars_area(self,RA=0,Dec=0,radius=10,**kwargs):
        """Get stars in a region of the sky

        Examples:
            # Get all stars around aldebaran in a radius of 5 degrees and with magnitudes between -1 and 4
            hyades = allstars.get_stars_area(RA=aldebaran.data.RA,Dec=aldebaran.data.Dec,
                                             radius=5,Mag=[-1,4])
        """
        stars = self.get_stars(
            RA=[RA-radius/15,RA+radius/15],
            Dec=[Dec-radius,Dec+radius],
            **kwargs
        )
        return stars
    
    def plot_stars(self,labels=False,pad=0,figargs=dict(),stargs=dict()):
        """Plot stars in a given area
        """
        plt.style.use('dark_background')

        # Figure
        dfigargs = dict(figsize=(5,5))
        dfigargs.update(figargs)
        fig,axs = plt.subplots(1,1,**dfigargs)

        # Axis
        axs.set_facecolor('black')

        # Scatter
        dstargs = dict(marker='*',color='y')
        dstargs.update(stargs)

        size_by_mag = Montu._linear_map([-1.5,5],[200,1])
        axs.scatter(15*self.data.RA,self.data.Dec,
                    s=size_by_mag(self.data.Mag),
                    **dstargs)
        
        # Labels
        if labels:
            for index in self.data.index:
                star = self.data.loc[index]
                star.fillna('',inplace=True)
                name = star.ProperName if star.ProperName != '' else star.BayerFlamsteed
                axs.text(15*star.RA,star.Dec,f'{name}',
                        color='w',fontsize=8)

        # Decoration
        axs.set_xlabel('RAJ2000 [deg]',fontsize=10)
        axs.set_ylabel('DecJ2000 [deg]',fontsize=10)
        
        # Range
        rang = max(15*((self.data.RA).max()-(self.data.RA).min()),
                (self.data.Dec).max()-(self.data.Dec).min())
        axs.margins(pad*rang)
        
        axs.grid(alpha=0.2)
        axs.axis('equal')
        fig.tight_layout()

        SET_PLT_DEFAULT_STYLE()
        return fig,axs


###############################################################
# Observing site
###############################################################
class ObservingSite(object):
    """Create an observing site

    Attributes:
        
        location: dictionary:
            Geodetic coordinates of the Observing site:
                lon: float [deg]: geodetic longitude
                lat: float [deg]: geodetic latitude
                elevation: float [km]: elevation

        conditions: dictionary:
            Physical conditions of the atmosphere on the observing site:
                pressure: surface pressure.
                temperature: surface temperature.
                relative_humidity: humidity of the air.
                obswl: observing band.

        planet: PlanetaryBody:
            Planetary body where the location is defined
    """
    def __init__(self,
                 planet=None,
                 lon=0,lat=0,height=0,
                 pressure=0,temperature=0,
                 relative_humidity=0,obswl=1):

        # Planet of the site
        self.planet = planet

        # Properties of the site
        self.lon = lon
        self.lat = lat
        self.elevation = height
        self.location = dict(lon=self.lon,lat=self.lat,elevation=self.elevation)

        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity
        self.obswl = obswl
        self.physical = dict(pressure=self.pressure*monu.hPa,temperature=self.temperature*monu.deg_C,
                               relative_humidity=self.relative_humidity,obswl=self.obswl*monu.micron),

        # Creation of site in astropy
        location_kwargs = dict(lon=lon*monu.deg,lat=lat*monu.deg,height=height*monu.km)
        self.astropy_site = astropy_EarthLocation(**location_kwargs)
        
        # Cartesian coordinates of the observer
        self.pos_fixedplanet = np.array(list(self.astropy_site.value))

        # Update site at J2000
        self.update_site(MonTime(0,format='tt',scale='tt'))

    def update_site(self,mtime):  
        """Update site
        """

        # Update epoch
        self.epoch = mtime

        # Rotation matrix from fixed frame on planet to J2000
        self.M_fixedplanet_J2000 = spy.pxform(self.planet.frameplanet,'J2000',self.epoch.tt)
        self.M_J2000_fixedplanet = spy.invert(self.M_fixedplanet_J2000)

        # Position of observer in J2000
        self.pos_J2000 = spy.mxv(self.M_fixedplanet_J2000,self.pos_fixedplanet)
        self.pos_Epoch = spy.mxv(self.epoch.M_J2000_Epoch,self.pos_J2000)

        # Update local reference frame vectors
        normal=spy.surfnm(self.planet.Re,self.planet.Re,self.planet.Rp,self.pos_fixedplanet)
        uy=spy.ucrss(np.array([0,0,1]),normal)
        ux=spy.ucrss(normal,uy)
        self.M_local_fixedplanet=np.array(np.vstack((ux,uy,normal)).transpose().tolist())
        self.M_fixedplanet_local=spy.invert(self.M_local_fixedplanet)

        # Position of site with respect to local system of reference
        self.pos_local = spy.mxv(self.M_fixedplanet_local,self.pos_fixedplanet)

        # Compute true sidereal time
        self.ltst = self.epoch.gtst + self.lon/15
        
###############################################################
# Planetary body
###############################################################
class PlanetaryBody(object):
    """Create a planetary body

    Examples: 
        earth = PlanetaryBody('399')
        mars = PlanetaryBody('mars')
        jupiter = PlanetaryBody('JUPITER')

    """

    @lru_cache()
    def __new__(cls,id):
        """This method is intended to avoid creating a new object with the same id
        Instead this method create a clone of the previously created object.
        """
        return super().__new__(cls)

    def __init__(self,body):

        # Names
        self.body = body.upper()
        if self.body in PLANETARY_IDS.keys():
            self.name = self.body.lower()
            self.id = str(PLANETARY_IDS[self.body])
        elif self.body in PLANETARY_NAMES.keys():
            self.id = str(self.body)
            self.name = PLANETARY_NAMES[self.id].lower()
        else:
            raise ValueError(f"Planet '{self.body}' not recognized, check variable PLANETARY_NAMES")
        self.capital = self.name[0].upper() + self.name[1:]

        # IAU frame of planet
        self.frameplanet = f'IAU_{self.name.upper()}'

        # Obtain object properties
        try:
            n,rs=spy.bodvrd(self.id,'RADII',3)
        except:
            n,rs=spy.bodvrd(self.id+'99','RADII',3)
        self.Re=rs[0]
        self.Rt=rs[1]
        self.Rp=rs[2]
        self.f=(self.Re-self.Rp)/self.Re

        # Create object from PyPlanet: try to avoid 
        exec(f'self.pymeeus_planet = pymeeus_{self.capital}')

        # Create Horizons object
        self.query_horizons = astropy_Horizons(id='4')

        # No predicted yet
        self.predict = False

    def calculate_sky_position(self,mtime=None,site=None,method='SPICE',
                               store=None,verbose=True):
        """Calculate position of a planet 

        Parameters:
            mtime: MonTime:
                Time when the position of the planet will be calculated.

            site: Site:
                Observing site w.r.t. position of the planet will be calculated.

            method:
                'Horizons': use astroquery.
                'SPICE': use SPICE and kernels.
                'VSOP87': use VSOP87 analytical theory.
                'all': all methods
        """
        self.data = {
            'datetime64':[mtime.obj_datejulian],
            'tt':[mtime.tt],
            'jtd':[mtime.jtd],
            'jed':[mtime.jed],

            }
        if store:
            if 'df' not in self.__dict__.keys():
                self.df = pd.DataFrame()
        
        if site is None:
            raise ValueError("No site selected")
        self.site = site

        # Update orientation of site
        Montu.vprint(verbose,f"Computing position of body '{self.name}' at epoch: jtd = {mtime.jtd} ")
        if mtime.tt != site.epoch.tt:
            Montu.vprint(verbose,f"Updating orientation of site (old time {site.epoch.datespice}, new time {mtime.datespice})")
            site.update_site(mtime)

        # Storing epoch
        self.epoch = mtime
            
        # Check if all methods are asked
        all_methods = True if method == 'all' else False

        # Check compute
        compute = False

        # Compute absolute coordinates (RA & Dec) at J2000 and epoch
        if method == 'Horizons' or all_methods:            

            # Using Horizons database
            Montu.vprint(verbose,"Method 'Horizons':")

            # Query ephemerides
            self.query_horizons.location = site.location
            self.query_horizons.epochs = self.epoch.jed # Julian Day as it was UTC
            self.ephemerides = self.query_horizons.ephemerides().to_pandas()
            ephemerides = self.ephemerides.loc[0]

            # Sky coordinates in J2000
            self.RAJ2000 = float(ephemerides.RA/15)
            self.DecJ2000 = float(ephemerides.DEC)

            # Distance
            self.distance = float(ephemerides.delta)
            self.sundistance = float(ephemerides.r)

            # Elongation and phase
            self.elongation = float(self.ephemerides.elong)
            self.phase = float(self.ephemerides.alpha_true)

            # Magnitude
            self.magnitude = self.pymeeus_planet.magnitude(self.sundistance,
                                                           self.distance,
                                                           self.phase*DEG)

            # Ecliptic coordinates in J2000
            uJ2000 = np.array([np.cos(self.DecJ2000*DEG)*np.cos(15*self.RAJ2000*DEG),
                               np.cos(self.DecJ2000*DEG)*np.sin(15*self.RAJ2000*DEG),
                               np.sin(self.DecJ2000*DEG)])
            ecJ2000 = spy.mxv(M_J2000_ECLIPJ2000,uJ2000)
            r,LonJ2000,LatJ2000 = spy.recrad(ecJ2000)
            self.LonJ2000 = LonJ2000*RAD
            self.LatJ2000 = LatJ2000*RAD

            # Sky coordinates @ Epoch
            self.RAEpoch = float(ephemerides.RA_app/15)
            self.DecEpoch = float(ephemerides.DEC_app)
            self.LonEpoch = float(ephemerides.ObsEclLon)
            self.LatEpoch = float(ephemerides.ObsEclLat)
            
            # Apparent coordinates @ site
            self.az = float(ephemerides.AZ)
            self.el = float(ephemerides.EL)
            
            # Compute hour angle
            self.HA = np.mod(np.arctan2(-np.sin(self.az*DEG)*np.cos(self.el*DEG)/np.cos(self.DecEpoch*DEG),
                                        (np.sin(self.el*DEG) - np.sin(self.DecEpoch*DEG)*np.sin(site.lat*DEG)) / \
                                        (np.cos(self.DecEpoch*DEG) * np.cos(site.lat*DEG)))*RAD/15,24)
            
            # Sidereal time
            self.tsa = np.mod(self.RAEpoch + self.HA,24)

            self._store_data('Horizons')
            self._show_position(verbose)

            compute = True

        if method == 'VSOP87' or all_methods:

            # Using VSOP87 semianalytical model implemented in PyEphem
            Montu.vprint(verbose,"Method 'VSOP87':")

            # Prepare object
            self.pyephem_planet = eval(f'pyephem.{self.capital}()')

            # Define observer
            self.pyephem_site = pyephem.Observer()
            self.pyephem_site.date = self.epoch.obj_pyephem
            self.pyephem_site.lat = f'{site.lat}'
            self.pyephem_site.lon = f'{site.lon}'
            self.pyephem_site.elevation = site.elevation
            self.pyephem_site.temp = site.temperature
            self.pyephem_site.pressure = site.pressure
            
            # Compute ephemerides
            self.pyephem_planet.compute(self.pyephem_site)
            self.RAEpoch = float(self.pyephem_planet.ra)*RAD/15
            self.DecEpoch = float(self.pyephem_planet.dec)*RAD

            # Distance
            self.distance = float(self.pyephem_planet.earth_distance)
            self.sundistance = float(self.pyephem_planet.sun_distance)
   
            # Elongation and phase
            self.elongation = abs(self.pyephem_planet.elong*RAD)
            illum = float(self.pyephem_planet.phase)
            self.phase = np.arccos(2*(illum/100)-1)*RAD
            self.magnitude = float(self.pyephem_planet.mag)

            # Ecliptic coordinates at Epoch
            uEpoch = np.array([np.cos(self.DecEpoch*DEG)*np.cos(15*self.RAEpoch*DEG),
                               np.cos(self.DecEpoch*DEG)*np.sin(15*self.RAEpoch*DEG),
                               np.sin(self.DecEpoch*DEG)])
            ecEpoch = spy.mxv(site.epoch.M_equatorial_ecliptic,uEpoch)
            r,LonEpoch,LatEpoch = spy.recrad(ecEpoch)
            self.LonEpoch = LonEpoch*RAD
            self.LatEpoch = LatEpoch*RAD
            
            # Compute elevation and azimuth
            self.el = self.pyephem_planet.alt*RAD
            self.az = self.pyephem_planet.az*RAD

            # Precess towards J2000
            self.RAJ2000,self.DecJ2000 = precession_equatorial(mtime.obj_pyplanet,pyplanets_Epoch(JED_2000),
                                                               pyplanets_Angle(15*self.RAEpoch),
                                                               pyplanets_Angle(self.DecEpoch))
            self.RAJ2000 = float(self.RAJ2000)/15
            self.DecJ2000 = float(self.DecJ2000)

            # Ecliptic coordinates in J2000
            uJ2000 = np.array([np.cos(self.DecJ2000*DEG)*np.cos(15*self.RAJ2000*DEG),
                               np.cos(self.DecJ2000*DEG)*np.sin(15*self.RAJ2000*DEG),
                               np.sin(self.DecJ2000*DEG)])
            ecJ2000 = spy.mxv(M_J2000_ECLIPJ2000,uJ2000)
            r,LonJ2000,LatJ2000 = spy.recrad(ecJ2000)
            self.LonJ2000 = LonJ2000*RAD
            self.LatJ2000 = LatJ2000*RAD
                        
            # Compute auxiliar coordinates
            self.HA = np.mod(np.arctan2(-np.sin(self.az*DEG)*np.cos(self.el*DEG)/np.cos(self.DecEpoch*DEG),
                                        (np.sin(self.el*DEG) - np.sin(self.DecEpoch*DEG)*np.sin(site.lat*DEG)) / \
                                        (np.cos(self.DecEpoch*DEG) * np.cos(site.lat*DEG)))*RAD/15,24)
            self.tsa = np.mod(self.RAEpoch + self.HA,24)

            self._store_data('VSOP87')
            self._show_position(verbose)

            compute = True

        if method == 'SPICE' or all_methods:

            # Using SPICE+pyplanets tools
            Montu.vprint(verbose,"Method 'SPICE':")

            # Retrieve positions in space
            site_planet_SSB_J2000,lt = spy.spkezr(site.planet.id,mtime.tt,'J2000','None','SSB')
            planet_SSB_J2000,lt = spy.spkezr(self.id,mtime.tt,'J2000','None','SSB')
            sun_SSB_J2000,lt = spy.spkezr('10',mtime.tt,'J2000','None','SSB')
            site_SSB_J2000 = site_planet_SSB_J2000[:3] + site.pos_J2000 

            # Celestial Coordinates at J2000
            planet_site_J2000 = planet_SSB_J2000[:3] - site_SSB_J2000
            r,RAJ2000,DECJ2000 = spy.recrad(planet_site_J2000)
            self.RAJ2000 = RAJ2000*RAD/15
            self.DecJ2000 = DECJ2000*RAD

            # Phase angle
            site_sun_J2000 = site_SSB_J2000 - sun_SSB_J2000[:3]
            planet_sun_J2000 = planet_SSB_J2000[:3] - sun_SSB_J2000[:3]

            self.elongation = np.arccos(-site_sun_J2000@planet_site_J2000/\
                ((site_sun_J2000@site_sun_J2000)**0.5*(planet_site_J2000@planet_site_J2000)**0.5))*RAD
            self.phase = np.arccos(planet_sun_J2000@planet_site_J2000/\
                ((planet_site_J2000@planet_site_J2000)**0.5*(planet_sun_J2000@planet_sun_J2000)**0.5))*RAD
            
            # Distance
            self.distance = (planet_site_J2000 @ planet_site_J2000)**0.5/AU
            self.sundistance = (planet_sun_J2000 @ planet_sun_J2000)**0.5/AU

            # Magnitude
            self.magnitude = self.pymeeus_planet.magnitude(self.sundistance,
                                                           self.distance,
                                                           self.phase*DEG)

            # Ecliptic coordinates J2000
            planet_site_EJ2000 = spy.mxv(M_J2000_ECLIPJ2000,planet_site_J2000)
            r,LonJ2000,LatJ2000 = spy.recrad(planet_site_EJ2000)
            self.LonJ2000 = LonJ2000*RAD
            self.LatJ2000 = LatJ2000*RAD

            # Celestial Coordinates at Epoch
            planet_site_Epoch = spy.mxv(self.epoch.M_J2000_Epoch,planet_site_J2000)
            r,RAplanet,DECplanet = spy.recrad(planet_site_Epoch)
            self.RAEpoch = RAplanet*RAD/15
            self.DecEpoch = DECplanet*RAD

            # Ecliptic coordinates at Epoch
            uEpoch = np.array([np.cos(self.DecEpoch*DEG)*np.cos(15*self.RAEpoch*DEG),
                               np.cos(self.DecEpoch*DEG)*np.sin(15*self.RAEpoch*DEG),
                               np.sin(self.DecEpoch*DEG)])
            ecEpoch = spy.mxv(site.epoch.M_equatorial_ecliptic,uEpoch)
            r,LonEpoch,LatEpoch = spy.recrad(ecEpoch)
            self.LonEpoch = LonEpoch*RAD
            self.LatEpoch = LatEpoch*RAD

            # Compute hour angle
            self.HA = np.mod(site.ltst - self.RAEpoch,24)
            
            # Compute elevation and azimuth
            self.el = np.arcsin(np.sin(self.DecEpoch*DEG)*np.sin(site.lat*DEG) + \
                                np.cos(self.DecEpoch*DEG)*np.cos(site.lat*DEG)*np.cos(self.HA*15*DEG))*RAD
            self.az = np.arctan2(-np.sin(self.HA*15*DEG)*np.cos(self.DecEpoch*DEG)/np.cos(self.el*DEG),
                                 (np.sin(self.DecEpoch*DEG) - np.sin(site.lat*DEG)*np.sin(self.el*DEG))/\
                                    (np.cos(site.lat*DEG)*np.cos(self.el*DEG)))*RAD
            self.az = np.mod(self.az,360)


            # Local sidereal time
            self.tsa = np.mod(self.RAEpoch + self.HA,24) 

            self._store_data('SPICE')
            self._show_position(verbose)

            compute = True

        if not compute:
            raise ValueError(f"Method '{method}' for computing ephemerides not recognized")

        if store:
            self.df=pd.concat([self.df,pd.DataFrame(self.data)])
            self.df.reset_index(drop=True,inplace=True)

    def _show_position(self,verbose):
        Montu.vprint(verbose,"\tPosition Epoch: ",
                     self.epoch.datestr,
                     self.epoch.jed)
        Montu.vprint(verbose,"\tCoordinates @ J2000: ")
        Montu.vprint(verbose,"\t\tEquatorial:",D2H(self.RAJ2000),D2H(self.DecJ2000))
        Montu.vprint(verbose,"\t\tEcliptic:",D2H(self.LonJ2000),D2H(self.LatJ2000))
        Montu.vprint(verbose,f"\tCoordinates @ Epoch : ")
        Montu.vprint(verbose,f"\t\tEquatorial:",D2H(self.RAEpoch),D2H(self.DecEpoch))
        Montu.vprint(verbose,f"\t\tEcliptic:",D2H(self.LonEpoch),D2H(self.LatEpoch))
        Montu.vprint(verbose,f"\tObserving conditions: ")
        Montu.vprint(verbose,f"\t\tDistance to site [au]: ",self.distance)
        Montu.vprint(verbose,f"\t\tDistance to sun [au]: ",self.sundistance)
        Montu.vprint(verbose,f"\t\tSolar elongation [deg]: ",D2H(self.elongation))
        Montu.vprint(verbose,f"\t\tPhase angle [deg]: ",D2H(self.phase))
        Montu.vprint(verbose,f"\t\tMagnitude: ",self.magnitude)
        Montu.vprint(verbose,f"\tOther properties: ")
        Montu.vprint(verbose,f"\t\tLocal true sidereal time: ",D2H(self.site.ltst))
        Montu.vprint(verbose,f"\t\tHour angle @ Epoch: ",D2H(self.HA))
        Montu.vprint(verbose,f"\t\tLocal coordinates @ Epoch: ",D2H(self.az),D2H(self.el))
        
    def ecliptic_longitude_advance(self,mtime,site,dt=1*HOUR,method='SPICE'):
        """Compute the rate of ecliptic longitude advance
        """

        # Time before
        self.calculate_sky_position(mtime-dt,site,method,verbose=0)
        EclLon_m_dt = self.LonEpoch

        # Time after 
        self.calculate_sky_position(mtime+dt,site,method,verbose=0)
        EclLon_p_dt = self.LonEpoch

        # Angle diff
        angle_diff = (EclLon_p_dt-EclLon_m_dt)

        # Angle differences: thanx ChatGPT!
        if angle_diff > 180:
            angle_diff -= 360
        elif angle_diff < -180:
            angle_diff += 360
        
        # Compute derivative using central difference algorithm
        dlondt = angle_diff/(2*dt)*DAY # Degrees per day

        return dlondt

    def reset_store(self):
        self.df = pd.DataFrame()

    def _store_data(self,metstr):
        self.data.update({
                # Generic
                f'RAJ2000':[self.RAJ2000],f'DecJ2000':[self.DecJ2000],
                f'RAEpoch':[self.RAEpoch],f'DecEpoch':[self.DecEpoch],
                f'LonJ2000':[self.LonJ2000],f'LatJ2000':[self.LatJ2000],
                f'LonEpoch':[self.LonEpoch],f'LatEpoch':[self.LatEpoch],
                f'tsa':[self.tsa],f'HA':[self.HA],f'az':[self.az],f'el':[self.el],
                f'site_distance':[self.distance],f'sun_distance':[self.sundistance],
                f'elongation':[self.elongation],f'phase':[self.phase],f'mag':[self.magnitude],
                # Method
                f'RAJ2000_{metstr}':[self.RAJ2000],f'DecJ2000_{metstr}':[self.DecJ2000],
                f'RAEpoch_{metstr}':[self.RAEpoch],f'DecEpoch_{metstr}':[self.DecEpoch],
                f'LonJ2000_{metstr}':[self.LonJ2000],f'LatJ2000_{metstr}':[self.LatJ2000],
                f'LonEpoch_{metstr}':[self.LonEpoch],f'LatEpoch_{metstr}':[self.LatEpoch],
                f'tsa_{metstr}':[self.tsa],f'HA_{metstr}':[self.HA],f'az_{metstr}':[self.az],f'el_{metstr}':[self.el],
                f'site_distance_{metstr}':[self.distance],f'sun_distance_{metstr}':[self.sundistance],
                f'elongation_{metstr}':[self.elongation],f'phase_{metstr}':[self.phase],f'mag_{metstr}':[self.magnitude],
            })

    def __str__(self):
        str = f"""Planetary Body
===============================
Names: 
    Basic name: {self.name}
"""
        return str

class Extra(object):
    def station_longitude_1(y):
            """Based on PyMeeus
            Valid only for [-2000,4000]
            """
            # Set some specific constants for Mars' opposition
            a = 2452097.382
            b = 779.936104
            m0 = 181.9573
            m1 = 48.705244
            k = round((365.2425 * y + 1721060.0 - a) / b)
            jde0 = a + k * b
            m = m0 + k * m1
            m = m*DEG
            t = (jde0 - 2451545.0) / 36525.0
            corr = (-37.079 + t * (-0.0009 + t * 0.00002) +
                    np.sin(m) * (-20.0651 + t * (0.0228 + t * 0.00004)) +
                    np.cos(m) * (14.5205 + t * (0.0504 - t * 0.00001)) +
                    np.sin(2.0 * m) * (1.1737 - t * 0.0169) +
                    np.cos(2.0 * m) * (-4.255 + t * (-0.0075 + t * 0.00008)) +
                    np.sin(3.0 * m) * (0.4897 + t * (0.0074 - t * 0.00001)) +
                    np.cos(3.0 * m) * (1.1151 + t * (-0.0021 - t * 0.00005)) +
                    np.sin(4.0 * m) * (-0.3636 + t * (-0.002 + t * 0.00001)) +
                    np.cos(4.0 * m) * (-0.1769 + t * (0.0028 + t * 0.00002)) +
                    np.sin(5.0 * m) * (0.1437 - t * 0.0004) +
                    np.cos(5.0 * m) * (-0.0383 - t * 0.0016))
            jde = jde0 + corr
            return jde

    def station_longitude_2(y):
            """Based on PyMeeus
            Valid only for [-2000,4000]
            """
            # Set some specific constants for Mars' opposition
            a = 2452097.382
            b = 779.936104
            m0 = 181.9573
            m1 = 48.705244
            k = round((365.2425 * y + 1721060.0 - a) / b)
            jde0 = a + k * b
            m = m0 + k * m1
            m = m*DEG
            t = (jde0 - 2451545.0) / 36525.0
            corr = (36.7191 + t * (0.0016 + t * 0.00003) +
                    np.sin(m) * (-12.6163 + t * (0.0417 - t * 0.00001)) +
                    np.cos(m) * (20.1218 + t * (0.0379 - t * 0.00006)) +
                    np.sin(2.0 * m) * (-1.636 - t * 0.019) +
                    np.cos(2.0 * m) * (-3.9657 + t * (0.0045 + t * 0.00007)) +
                    np.sin(3.0 * m) * (1.1546 + t * (0.0029 - t * 0.00003)) +
                    np.cos(3.0 * m) * (0.2888 + t * (-0.0073 - t * 0.00002)) +
                    np.sin(4.0 * m) * (-0.3128 + t * (0.0017 + t * 0.00002)) +
                    np.cos(4.0 * m) * (0.2513 + t * (0.0026 - t * 0.00002)) +
                    np.sin(5.0 * m) * (-0.0021 - t * 0.0016) +
                    np.cos(5.0 * m) * (-0.1497 - t * 0.0006))
            jde = jde0 + corr
            return jde   

###############################################################
# Individual modules
###############################################################
"""
from montu.time import *
from montu.planets import *
from montu.stars import *
"""

###############################################################
# Tests 
###############################################################
#"""
if __name__ == '__main__':
    print('Montu Python Test Suite')

    # Util tools
    print("Testing vprint:")
    Montu.vprint(False,"Hidden message")
    Montu.vprint(True,"Visible message")

    # File path
    print("Testing data path:")
    print(Montu._data_path('frame.tk'))

    # Testing files
    print("Testing kernels")
    Montu._wget('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls',
                '/tmp/leaps.tls',verbose=True)
    
    Montu.load_kernels({'pck00011.tpc':''})
    Montu.load_kernels()
    Montu.load_kernels()

    # Testing stars
    allstars = Stars()
    print("Number of stars in database: ",len(allstars.data))

    # Get a single stars
    aldebaran = allstars.get_stars(ProperName='Aldebaran')
    print("Aldebaran J2000: ",aldebaran.data)

    # All visible stars in the sky
    visible = allstars.get_stars(Mag=[-2,6.5])
    print("Number of visible stars in the sky",visible.number)

    # All visible stars with declination less than 1 deg
    equator = allstars.get_stars(Mag=[-2,6.5],Dec=[-1,1])
    print("Number of visible stars near the equator",equator.number)

    # Get Hyades
    hyades = allstars.get_stars_area(RA=aldebaran.data.RA,Dec=aldebaran.data.Dec,
                                     radius=5,Mag=[-1,4])
    print("Number of Hyades:",hyades.number)

    # Plot Hyades
    hyades = allstars.get_stars_area(RA=aldebaran.data.RA,Dec=aldebaran.data.Dec,
                                     radius=10,Mag=[-1,10])
    #fig = hyades.plot_stars(pad=0.0,labels=False,figargs=dict(figsize=(8,8)))
    #plt.show()

    # Dates
    print("Testing dates")
    print("Reference date: J2000")
    mtime = MonTime('2000-01-01 12:00:00.00',format='iso',scale='tt',proleptic=True,verbose=True)
    print()
    mtime = MonTime('2000-01-01 12:00:00.00',format='iso',scale='utc',proleptic=True,verbose=True)
    print()
    mtime = MonTime('2000-01-01 12:00:00.00',format='iso',scale='tt',proleptic=False,verbose=True)
    print()
    mtime = MonTime('2000-01-01 12:00:00.00',format='iso',scale='utc',proleptic=False,verbose=True)

    print("Ancient date: -2500")
    mtime = MonTime('-2500-01-01 12:00:00.00',format='iso',scale='utc',proleptic=True,verbose=True)
    print()
    mtime = MonTime('-2500-01-01 12:00:00.00',format='iso',scale='utc',proleptic=False,verbose=True)
    print()
    mtime = MonTime('-2500-01-01 12:00:00.00',format='iso',scale='tt',proleptic=True,verbose=True)
    print()
    mtime = MonTime('-2500-01-01 12:00:00.00',format='iso',scale='tt',proleptic=False,verbose=True)
    print()

    print("Initialize with terrestrial time")
    mtime = MonTime(807954,format='jd',scale='tt',proleptic=True,verbose=True)
    print()
    mtime = MonTime(807954,format='jd',scale='utc',proleptic=True,verbose=True)
    
#"""

