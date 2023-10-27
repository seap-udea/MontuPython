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
import regex
import math
import inspect

# Basic packages
import numpy as np
import spiceypy as spy
import pandas as pd

# Class utilities
from functools import lru_cache
from tabulate import tabulate

# Avoid warnings
import warnings
warnings.filterwarnings("ignore")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Atrotools
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Pyephem
import ephem as pyephem

# PyPlanets and pyMeeus
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

# Astropy
from astropy.time import Time as astropy_Time
from astroquery.jplhorizons import Horizons as astropy_Horizons
from astropy.coordinates import EarthLocation as astropy_EarthLocation

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Numerical toolstrotools
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Scipy
from scipy.interpolate import interp1d

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# REMOVE THIS PACKAGES WHEN DONE
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
from datetime import datetime
import matplotlib.pyplot as plt

# Astropy
from astropy import constants as astropy_constants
from astropy import units as monu # Montu Units are astropy units
from astropy.coordinates import AltAz as astropy_AltAz
from astropy.coordinates import SkyCoord as astropy_SkyCoord
from astropy.coordinates import Angle as astropy_Angle

###############################################################
# Global settings of the package
###############################################################
# Set precision
np.set_printoptions(precision=17)
pd.set_option("display.precision",17)

###############################################################
# Constants
###############################################################
# Numerical Constants
RAD = 180/np.pi
DEG = 1/RAD
UARCSEC = 1e-6/3600 # mircroarsec in degrees
MARCSEC = 1e-3/3600 # milliarcsec in degrees 

# Phsyical
AU = 149597870.700 # km, value in Horizons
CSPEED = 299792.458 # km/s, value in Horizons

# Astronomical
MIN = 60 # s
HOUR = 60*MIN # s
DAY = 86400 # s
CALYEAR = 365*DAY # s, calendar year
JULYEAR = 365.25*DAY # s, julian year
YEAR = 365.25*DAY # s
CENTURY = 100*YEAR # s
MILLENIUM = 10*YEAR # s
TAI_TO_SID = 1.00273781191135448 # Sidereal seconds / Uniform seconds
LETTERS=dict(A=1,B=2,C=3,D=4,E=5,F=6)

# Required kernels
"""This dictionaries describe the kernels the package require to compute planetary positions

If the dictionary is blank it means that the kernel is in the data directory.

The leapseconds kernes should be updated periodically.
"""
BASIC_KERNELS = {
    'naif0012.tls':'',
    'frame.tk':'',
    'pck00011.tpc':'',
    'earth_assoc_itrf93.tf':''
}
PRECISION_KERNELS = {
    'latest_leapseconds.tls':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls',
    'de441_part-1.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp',
    'de441_part-2.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp',
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

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# REMOVE THIS CONSTANTS WHEN DONE
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LEGACY = True
MAIN = False

SIDEREAL_YEAR = 365.256363004*DAY # d, J2000
TROPICAL_YEAR = 365.242190402*DAY # d, J2000

# Rotation matrix between J2000 and ECLIPJ2000
M_J2000_ECLIPJ2000 = spy.pxform('J2000','ECLIPJ2000',0)

# Historical
PYEPHEM_JD_REF = 2415020.0
PYEPHEM_MJD_2000 = 36525.0
JED_2000 = 2451545.0

# JD PRECISION
"""
Rounding errors may produce strange figures when dealing with 
Julian Days. We round-up all JD values to JD_PRECISION_FIGURES
to avoid this artifacts. 7 figures correspond to 0.01 seconds.
"""
JD_PRECISION_FIGURES = 8

# Correction for JED
jed_correction_data = np.loadtxt(Montu._data_path('corrections_dt.dat'))
JED_CORRECTION = interp1d(jed_correction_data[:,0],jed_correction_data[:,1])

###############################################################
# Montu Python Util Class
###############################################################
class Seba(object):
    """This is the general class for all celestial objects in MontuPython

    The name of the class come from the word for star in ancient egyptian, 
    seba /sbꜣ/.
    """
    def __init__(self):
        pass

class Sebau(object):
    """This is the general class for a set of celestial objects in MontuPython

    The name of the class is the plural for star in ancient egyptian, 
    sebau /sbꜣw/.

    Basic attributes:
        number: number of objects
    """
    def __init__(self):
        pass

    def __repr__(self):
        return str(self.data.head(5))

    def __str__(self):
        str = self.data.head(5)
        return str

class Dictobj(object):
    """Convert a dictionary to an object

    Examples:
        ob = Dictobj(a=2,b=3)
        print(ob.a,ob.b)
        ob = Dictobj(dict=dict(a=2,b=3))
        print(ob.a,ob.b)
        ob = Dictobj(dict={'a':2,'b':3})
        print(ob.a,ob.b)
    """

    def __init__(self, **kwargs):
        if 'dict' in kwargs.keys():
            kwargs.update(kwargs['dict'])
        for key, value in kwargs.items():
            if key == 'dict':continue
            setattr(self, key, value)
        
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
# Montu Python Util Class
###############################################################
class Util(object):

    def vprint(verbose,*args):
        """Print messages in verbose mode
        """
        if verbose:
            print(*args)

    def arange(start, stop, step=1, endpoint=True):
        """Same as np.arange but including endpoint

        Source: https://stackoverflow.com/a/68551927
        """
        arr = np.arange(start, stop, step)
        if endpoint and arr[-1]+step==stop:
            arr = np.concatenate([arr,[stop]])
        return arr

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

    def load_kernels(kernels=BASIC_KERNELS,dir='montmp/',verbose=False):
        
        # Check if dir exists
        if not os.path.exists(dir):
            os.system(f"mkdir -p {dir}")

        # Load kernel
        for kernel,item in kernels.items():

            #
            if kernel in KERNELS_LOADED.keys():
                Util.vprint(verbose,f"Kernel {kernel} already loaded, skipping")
                continue

            # Local kernel
            if len(item) == 0:
                if verbose:print(f"Loading local kernel {kernel}")
                kernel_file = Util._data_path(kernel)
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
                Util._wget(item,kernel_path)
            
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

    def table_df(df,format='github'):
        """Present a DataFrame in a tabular form

        format: string, default = 'github':
            Format of the table.

            Other formats: “plain”,“simple”,“github”,“grid”,“fancy_grid”,“pipe”,
            “orgtbl”,“jira”,“presto”,“pretty”,“psql”,“rst”,“mediawiki”,“moinmoin”,
            “youtrack”,“html”,“latex”,“latex_raw”,“latex_booktabs”,“textile”, 
        """
        print(tabulate(df,headers='keys',tablefmt=format))

    def dec2hex(dec,string=True):

        dec = float(dec)
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

    def get_methods(my_class):
        """Get a list of the methods for class my_class
        """
        return sorted([member[0] for member in inspect.getmembers(my_class) if '__' not in member[0]])

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
        Util.vprint(verbose,f"Downloading {filename} from {url} [size = {total_size}]")

        # Initialize the progress bar
        progress_bar = tqdm.tqdm(total=total_size, unit='B', unit_scale=True)

        with open(filename, 'wb') as file:
            for data in response.iter_content(chunk_size=1024):
                file.write(data)
                progress_bar.update(len(data))

        progress_bar.close()

    def _linear_map(mapped,observed):
        a = (observed[1]-observed[0])/(mapped[1]-mapped[0])
        b = observed[0] - a*mapped[0]
        map = lambda x:a*x+b
        return map

# Aliases
D2H = Util.dec2hex
PRINTDF = Util.print_df
TABLEDF = Util.table_df

###############################################################
# Main class
###############################################################
class Time(object):
    """Create a time-frame object
    
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

        proleptic: boolean, default = True

            Proleptic gregorian correspond to the case when the Gregorian calendar
            is extended before the adoption date at 1582-10-15. When proleptic = False, the Julian 
            calendar is used before the adoption date.

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


            obj_astrotyme: Time:
                Date in astropy format.
            
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
                 date=None,
                 format='iso',
                 scale='utc',
                 calendar='proleptic'):
        

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

        # Set calendar assumed
        self.calendar = calendar
        
        if format == 'iso':
            # Date was provided as a string

            if calendar=='proleptic':

                # Parse string
                self._parse_datestr(date)

                # Calculated deltat
                deltat = pymeeus_Epoch.tt2ut(self.components[0]*self.components[1],self.components[2])

                # Convert to terrestrial time as if date was given in TT scale
                et = spy.utc2et(self.datespice)
                deltat_leaps = spy.deltet(et,'ET')
                et -= deltat_leaps
                
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
                jed = round(spy.unitim(et,'ET','JED'),JD_PRECISION_FIGURES)
                jtd = round(spy.unitim(tt,'ET','JED'),JD_PRECISION_FIGURES)

                """
                There is an error of between 0.5 and 10 seconds for years 
                between 300 c.e. and 1582 c.e. of unknown origin and it seems
                to come from the SPICE algorithms or the algorithms calculating
                deltat in PyMeeus. This code correct this effect.
                """
                year = self.components[0]*self.components[1]
                if 300<=year<=1582:
                    jed_correction = JED_CORRECTION(year)
                    tt -= jed_correction

            elif calendar=='mixed':

                # Parse string
                self._parse_datestr(date)

                # Calculated deltat
                deltat = pymeeus_Epoch.tt2ut(self.components[0]*self.components[1],self.components[2])

                # Convert from components to pymeeus epoch which receives date in mixed calendar
                args = (self.components[0]*self.components[1],
                        self.components[2],
                        self.components[3]+\
                        self.components[4]/24+\
                        self.components[5]/(24*60)+\
                        self.components[6]/(24*60*60)+\
                        self.components[7]/(24*60*60*1e6))
                pymeeus_epoch = pymeeus_Epoch(*args)
                jd = pymeeus_epoch.jde()
                et = spy.unitim(jd,'JED','ET')
                
                # According to scale choose terrestrial time
                if scale == 'tt':
                    tt = et
                    jtd = round(jd,JD_PRECISION_FIGURES)
                    et = tt - deltat
                    jed = round(jtd - deltat/DAY,JD_PRECISION_FIGURES)
                else:
                    et = et
                    jed = round(jd,JD_PRECISION_FIGURES)
                    tt = et + deltat
                    jtd = round(jed + deltat/DAY,JD_PRECISION_FIGURES)

            else:
                raise ValueError("Calendar '{calendar}' not recognzed. Use 'proleptic' or 'mixed'.")

            # Initialize object according to tt
            self.update_time(tt,format='tt',scale='tt')

        else:
            # Initialize object according to date and format and 
            self.update_time(date,format,scale)

    def _parse_datestr(self,date):
        """Parse date string
        """

        # Default format
        style = 'astro' # Default style of input string 
        self.datepro = date # Default date

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
            self.datepro = re.sub('bce[a-z\s]*(\d+)',subs1,self.datepro.lower())
            self.datepro = re.sub('(\d+)\s*b[\.]*c[\.]*e[\.]*\s*',subs2,self.datepro.lower())

        # Create calendar and datetime object
        self.obj_datetime64 = np.datetime64(self.datepro)
        self.components = Util.dt2cal(np.datetime64(self.datepro.strip('-')),
                                       bce=self.bce)
        
        # Generate info
        self.year = self.components[0]*self.components[1]
        self.month = self.components[2]
        self.day = self.components[3]
        self.hour = self.components[4]
        self.minute = self.components[5]
        self.second = self.components[6]
        self.usecond = self.components[7]
        
        # Adjust SPICE string according to epoch
        if self.bce:
            self.datespice = f'{-self.year+1:04d} B.C. {self.month:02d}-{self.day:02d} {self.hour:02d}:{self.minute:02d}:{self.second:02d}.{self.usecond:02d}'
        elif 0<self.year<1000:
            self.datespice = f'{self.year:04d} A.D. {self.month:02d}-{self.day:02d} {self.hour:02d}:{self.minute:02d}:{self.second:02d}.{self.usecond:02d}'
        else:
            self.datespice = self.datepro
        
    def update_time(self,time=None,format='tt',scale='tt'):
        """Update time object according to terrestrial time
        
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
            jd = round(spy.unitim(et,'ET','JED'),JD_PRECISION_FIGURES)
        else:
            raise AssertionError(f"Format '{format}' not recognized (valid 'iso', 'tt', 'jd')")

        # Initialize Epoch
        pymeeus_epoch = pymeeus_Epoch(jd)
        year,month,day = pymeeus_epoch.get_date()
        self.isjulian = pymeeus_Epoch.is_julian(year,month,day)
        self.deltat = pymeeus_Epoch.tt2ut(year,month)
        self.bce = True if year<=0 else False
        
        # Terrestrial time
        et = spy.unitim(jd,'JED','ET')
        if scale == 'tt':
            self.jtd = round(jd,JD_PRECISION_FIGURES)
            self.tt = et
            self.jed = round(jd - self.deltat/DAY,JD_PRECISION_FIGURES)
            self.et = self.tt - self.deltat
        else:
            self.jed = round(jd,JD_PRECISION_FIGURES)
            self.et = et
            self.jtd = round(self.jed + self.deltat/DAY,JD_PRECISION_FIGURES)
            self.tt = self.et + self.deltat

        # Create pyplanet epoch: you need to provide jd with no deltat correction: is internal
        self.obj_pyplanet = pyplanets_Epoch(self.jed)

        # Create astrotime: you need the tdb time
        self.obj_astrotime = astropy_Time(self.jtd,format='jd',scale='tdb')
        
        # PyEphem Date: you need to provide jd with no deltat correction: is internal
        self.obj_pyephem = pyephem.Date(self.jed - PYEPHEM_JD_REF)
        
        # String for datemixed
        pyephem_str = f'{self.obj_pyephem}'.strip('-')
        parts = pyephem_str.split(' ')
        cals = [int(p) for p in parts[0].split('/')] +[int(p) for p in parts[1].split(':')] 
        if self.bce:
            # Adjust year if bce
            cals[0] -= 1
            cals[0] *= -1
        self.datemix = f'{cals[0]}-{cals[1]:02d}-{cals[2]:02d} {cals[3]:02d}:{cals[4]:02d}:{cals[4]:02d}'
        
        # Replace month name
        MONTH_ABREVS = dict(JAN=1,FEB=2,MAR=3,APR=4,MAY=5,JUN=6,JUL=7,AUG=8,SEP=9,OCT=10,NOV=11,DEC=12)

        # Set string from terrestrial time
        self.datespice = spy.et2utc(self.et+spy.deltet(self.et,'ET'),'C',4)
        datestr = self.datespice

        # Converting from 
        sub_bc = lambda m:f'-{int(m.group(1))-1:04d}-{MONTH_ABREVS[m.group(2)]:02d}-'
        sub_ad = lambda m:f'{int(m.group(1)):04d}-{MONTH_ABREVS[m.group(2)]:02d}-'
        sub_nm = lambda m:f'-{MONTH_ABREVS[m.group(1)]:02d}-'
        if self.bce:
            datestr = re.sub('(\d+)\s*B.C.\s*(\w+)\s*',sub_bc,datestr)
        elif 'A.D.' in datestr:
            datestr = re.sub('(\d+)\s*A.D.\s*(\w+)\s*',sub_ad,datestr)
        else:
            datestr = re.sub('\s+(\w+)\s+',sub_nm,datestr)
        
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
        return Util.string_difference(self._get_signature(),mtime._get_signature())

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
Date in proleptic UTC: {self.datepro}
Date in mixed UTC: {self.datemix}
Date in SPICE format: {self.datespice}
General:
    Components: {self.components}
    Is bce: {self.bce}
    Is Julian: {self.isjulian}
Uniform scales:
    Terrestrial time:
        tt: {self.tt}
        jtd: {self.jtd}
    UTC time:
        et: {self.et}
        jed: {self.jed}
    Delta-t = TT - UTC = {self.deltat}
Objects:
    Date in datetime64 format: {self.obj_datetime64}
    Date in PyPlanet Epoch: {self.obj_pyplanet}
    Date in PyEphem Epoch: {self.obj_pyephem}
    Date in AstroPy Time: {self.obj_astrotime}
Astronomical properties at Epoch:
    True obliquity of ecliptic: {D2H(self.epsilon,1)}
    True nutation longitude: {D2H(self.delta_psi,1)}
    Greenwhich Meridian Sidereal Time: {D2H(self.gtst,1)}
"""
        return str
    
    def __repr__(self) -> str:
        str = f"Time('{self.datepro}'/'{self.datemix}')"
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
    
    @staticmethod
    def tai_to_sid(mtime):
        """Compute the instantaneous rate of sidereal days to tai days at epoch mtime
        """
        # Present
        epsilon = float(pyplanets_true_obliquity(mtime.obj_pyplanet))
        delta_psi = float(pyplanets_nutation_longitude(mtime.obj_pyplanet))
        gtst1 = 24*mtime.obj_pyplanet.apparent_sidereal_time(epsilon,delta_psi)

        # Future
        mtime2 = mtime + 1*DAY
        epsilon = float(pyplanets_true_obliquity(mtime2.obj_pyplanet))
        delta_psi = float(pyplanets_nutation_longitude(mtime2.obj_pyplanet))
        gtst2 = 24*mtime2.obj_pyplanet.apparent_sidereal_time(epsilon,delta_psi)

        # Advance
        tai2sid = 1 + (gtst2 - gtst1)/24
        return tai2sid
    
    @staticmethod
    def set_time_ticks(ax):
        """Set xticks as Time objects
        """
        tts = ax.get_xticks()
        xlabels = []
        for tt in tts:
            mtime = Time(tt)
            xlabels += [f'{mtime.year}']
        ax.set_xticklabels(xlabels)

STELLAR_CATALOGUE = 'montu_stellar_catalogue_v37.csv'

"""References
    Styles: https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
    #sphx-glr-gallery-style-sheets-style-sheets-reference-py
"""
PLT_DEFAULT_STYLE = 'default' # others: ggplot, default, classic
SET_PLT_DEFAULT_STYLE = lambda:plt.style.use(PLT_DEFAULT_STYLE)

###############################################################
# Stars Class
###############################################################
class Stars(Sebau):
    """Stellar catalogue

    Initialization parameters:

        data: pandas.Dataframe, default = None:
            Pandas dataframe containing the stars.

        filename: string, default = None:
            File containing the database with stars.
            If None it uses the official MontuCatalogue, namely 
            montu_stellar_catalogue_x.csv (where x is the version of the catalogue)

    Attributes:
        
        data: pandas.DataFrame:
            Data containing the information on stars.

        number: int:
            Number of stars in data.

    Methods:
        get_stars:
            Get a subset of the stars in catalogue.

        get_stars_around:
            Get stars around a given position in the sky.

            
    """

    def __init__(self,data=None,filename=None):

        if data is not None:
            # Load data for stars from a dataframe already loaded
            self.data = copy.deepcopy(data)
            
        elif filename:
            # Load data from a file
            self.data = pd.read_csv(filename)

        else:
            # Load data from the database provided with package
            self.data = pd.read_csv(Util._data_path(STELLAR_CATALOGUE,check=True))

        self.number = len(self.data)

    def get_stars(self,**args):
        """Filter stars by criteria

        Examples:
            # Get a single stars
            aldebaran = allstars.get_stars(ProperName='Aldebaran')

            # All visible stars in the sky
            visible = allstars.get_stars(Mag=[-2,6.5])

            # All visible stars with declination less than 1 deg in absolute value
            equator = allstars.get_stars(Mag=[-3,6.5],Dec=[-1,1])
        """

        # If no args get all stars in data base
        if len(args)==0:
            return self
        
        # If args provided it will try to filter database according to conditions
        cond = np.array([True]*len(self.data))
        for key,item in args.items():
            if key == 'suffix':continue
            if isinstance(item,list):
                min = float(item[0])
                max = float(item[1])
                cond = (self.data[key]>=min)&(self.data[key]<=max)&(cond)
            elif isinstance(item,tuple):
                cond_or = np.array([False]*len(self.data))
                for it in item:
                    cond_or = (self.data[key]==it)|cond_or
                cond = (cond_or)&(cond)
            else:
                cond = (self.data[key]==item)&(cond)
    
        return Stars(self.data[cond])
    
    def get_stars_around(self,
                         center=[0,0],radius=10,
                         coords=['RAJ2000','DecJ2000'],**kwargs):
        """Get stars around a point in the sky

        Parameters:
            center: list, default = [0,0] 
                Center around which the region will be extracted.

            radius: float, default = 10 [same units as center]:
                Radius of the region

            coords: list of strings, default = ['RAJ2000','DecJ2000'];
                Name of coordinates on which the center is calculated.

        Return: 
            stars: Stars:
                Stars in area.

        Examples:
            # Get all stars around aldebaran in a radius of 5 degrees and with magnitudes between -1 and 4
            hyades = stars.get_stars_around(center=[aldebaran.data.RAJ2000,aldebaran.data.DecJ2000],radius=15,Vmag=[-1,4])
        """
        kwargs.update({
            coords[0]:[float(center[0]-radius/15),float(center[0]+radius/15)],
            coords[1]:[float(center[1]-radius),float(center[1]+radius)],
        })
        stars = self.get_stars(**kwargs)
        return stars
    
    def plot_stars(self,coords=['RAJ2000','DecJ2000'],
                   labels=True,pad=0,figargs=dict(),stargs=dict()):
        """Plot all stars in data.

        Parameters:
            coords: list of strings, default = ['RAJ2000','DecJ2000']
                Coordinates used in representation.

            labels: Boolean, default = True:
                Do you want to see stellar labels.

            figargs: dictionary:
                Additional options for the figure.
            
            starargs: dictionary:
                Additional options for the scatter.

        Return:
            fig,axs: 
                Figure and axes.
        """
        # Black background
        plt.style.use('dark_background')

        # Create figure
        dfigargs = dict(figsize=(8,8))
        dfigargs.update(figargs)
        fig,axs = plt.subplots(1,1,**dfigargs)

        # Axis
        axs.set_facecolor('black')

        # Scatter
        dstargs = dict(marker='*',color='y')
        dstargs.update(stargs)

        size_by_mag = Util._linear_map([-1.5,5],[200,1])
        axs.scatter(15*self.data[coords[0]],
                    self.data[coords[1]],
                    s=size_by_mag(self.data.Vmag),
                    **dstargs)
        
        # Labels
        fontsize = Util._linear_map([6,-2],[4,14])
        if labels:
            for index in self.data.index:
                star = self.data.loc[index]
                star.fillna('',inplace=True)
                name = star.Name
                """
                axs.text(star[coords[0]],star[coords[1]],f'{name}',
                         color='w',fontsize=fontsize(star.Vmag))
                         """
                axs.annotate(f'{name}',xy=[15*star[coords[0]],star[coords[1]]],xycoords='data',
                             xytext=[5,5],textcoords='offset points',
                             fontsize=fontsize(star.Vmag))

        # Decoration
        axs.set_xlabel(f'{coords[0]} [hh:mm]',fontsize=10)
        axs.set_ylabel(f'{coords[1]} [deg]',fontsize=10)
        
        # Range
        rang = max(((15*self.data[coords[0]]).max()-(15*self.data[coords[0]]).min()),
                   (self.data[coords[1]]).max()-(self.data[coords[1]]).min())
        axs.margins(pad*rang)
        
        axs.grid(alpha=0.2)
        axs.axis('equal')
        fig.tight_layout()

        # Change tick labels
        ra_ticks = axs.get_xticks()
        ra_tick_labels = []
        for ra in ra_ticks:
            comps = D2H(ra/15,string=False)
            ra_tick_labels += [f'{int(comps[0]):02d}:{comps[1]:02d}']
        axs.set_xticklabels(ra_tick_labels)

        dec_ticks = axs.get_yticks()
        dec_tick_labels = []
        for dec in dec_ticks:
            comps = D2H(dec,string=False)
            dec_tick_labels += [f'{int(comps[0]):02d}:{comps[1]:02d}']
        axs.set_yticklabels(dec_tick_labels,rotation=90)

        # Montu water mark
        Util.montu_mark(axs)

        SET_PLT_DEFAULT_STYLE()
        return fig,axs
    
    def __repr__(self):
        repr = f"Stars(number={len(self.data)})"
        return repr

    def __str__(self):
        desc = str(self.data)
        return desc

###############################################################
# Planetary body
###############################################################
class Planet(Seba):
    """Create a planetary body

    Attributes:
        body ('MARS'), name ('mars'), id ('4'), capital('Mars'): string

        frameplanet: string:
            Name of reference frame fixed on planet

    Examples: 
        earth = Planet('399')
        mars = Planet('mars')
        jupiter = Planet('JUPITER')

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
        try:
            exec(f'self.pymeeus_planet = pymeeus_{self.capital}')
        except:
            self.pymeeus_planet = None

        # Create Horizons object
        self.query_horizons = astropy_Horizons(id='4')

        # No predicted yet
        self.predict = False

        # Data: if None it indicates that no where_among_stars with store = True has been ran
        self.data = None

    def where_among_stars(self,mtime=None,site=None,
                          store=None,propermotion=True,
                          verbose=False):
        """Calculate position of a planet 

        Parameters:
            mtime: Time:
                Time when the position of the planet will be calculated.

            site: Site:
                Observing site w.r.t. position of the planet will be calculated.
        """
        self.position = {
            'jed':[mtime.jed],
            'tt':[mtime.tt],
            }
        if store:
            if 'data' not in self.__dict__.keys():
                self.data = pd.DataFrame()
        
        if site is None:
            raise ValueError("No site selected")
        self.site = site

        # Calculate coordinates
        position = self._calculate_coordinates(mtime.tt,site)
        self.RAJ2000,self.DecJ2000,self.LonJ2000,self.LatJ2000 = list(position)

        # Phase angle
        site_sun_J2000 = self.site_SSB_J2000 - self.sun_SSB_J2000[:3]
        planet_sun_J2000 = self.planet_SSB_J2000[:3] - self.sun_SSB_J2000[:3]

        self.elongation = np.arccos(-site_sun_J2000@self.planet_site_J2000/\
            ((site_sun_J2000@site_sun_J2000)**0.5*(self.planet_site_J2000@self.planet_site_J2000)**0.5))*RAD
        self.phase = np.arccos(planet_sun_J2000@self.planet_site_J2000/\
            ((self.planet_site_J2000@self.planet_site_J2000)**0.5*(planet_sun_J2000@planet_sun_J2000)**0.5))*RAD
        
        # Distance
        self.distance = (self.planet_site_J2000 @ self.planet_site_J2000)**0.5/AU
        self.sundistance = (planet_sun_J2000 @ planet_sun_J2000)**0.5/AU

        # Magnitude
        try:
            self.magnitude = self.pymeeus_planet.magnitude(self.sundistance,
                                                            self.distance,
                                                            self.phase*DEG)
        except:
            self.magnitude = 33
        
        # Calculate proper motion
        if propermotion == True:
            dpdt,d2pdt2 = self._calculate_proper_motion(mtime.tt,site,position=position)
            self.pmRA,self.pmDec,self.pmLon,self.pmLat = list(dpdt)
            self.paRA,self.paDec,self.paLon,self.paLat = list(d2pdt2)
        else:
            self.pmRA,self.pmDec,self.pmLon,self.pmLat = [0]*4
            self.paRA,self.paDec,self.paLon,self.paLat = [0]*4

        self._store_sky_position()
        if store:
            self.data=pd.concat([self.data,self.position])
            self.data.reset_index(drop=True,inplace=True)

    def _store_sky_position(self):
        self.position.update({
                f'RAJ2000':[self.RAJ2000],f'DecJ2000':[self.DecJ2000],
                f'pmRA':[self.pmRA],f'paRA':[self.paRA],f'pmDec':[self.pmDec],f'paDec':[self.paDec],
                f'LonJ2000':[self.LonJ2000],f'LatJ2000':[self.LatJ2000],
                f'pmLon':[self.pmLon],f'paLon':[self.paLon],f'pmLat':[self.pmLat],f'paLat':[self.paLat],
                f'site_distance':[self.distance],f'sun_distance':[self.sundistance],
                f'elongation':[self.elongation],f'phase':[self.phase],f'mag':[self.magnitude],
                f'XJ2000':[self.planet_site_J2000[0]],f'YJ2000':[self.planet_site_J2000[1]],f'ZJ2000':[self.planet_site_J2000[2]],
                f'VXJ2000':[self.vplanet_site_J2000[0]],f'VYJ2000':[self.vplanet_site_J2000[1]],f'VZJ2000':[self.vplanet_site_J2000[2]],
            })
        self.position = pd.DataFrame(self.position)

    def _compare_positions(self,mtime=None,site=None,
                     method='SPICE',store=None,
                     verbose=False):
        """Calculate position of a planet 

        Parameters:
            mtime: Time:
                Time when the position of the planet will be calculated.

            site: Site:
                Observing site w.r.t. position of the planet will be calculated.

            method:
                'Horizons': use astroquery.
                'SPICE': use SPICE and kernels.
                'VSOP87': use VSOP87 analytical theory.
                'all': all methods
        """
        self.position = {
            'datepro':[mtime.datepro],
            'datemix':[mtime.datemix],
            'datetime64':[mtime.obj_datetime64],
            'tt':[mtime.tt],
            'jtd':[mtime.jtd],
            'jed':[mtime.jed],
            }
        if store:
            if 'data' not in self.__dict__.keys():
                self.data = pd.DataFrame()
        
        if site is None:
            raise ValueError("No site selected")
        self.site = site

        Util.vprint(verbose,f"Computing position of body '{self.name}' at epoch: jtd = {mtime.jtd} ")
        
        # Update orientation of site
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
            Util.vprint(verbose,"Method 'Horizons':")

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
            Util.vprint(verbose,"Method 'VSOP87':")

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
            Util.vprint(verbose,"Method 'SPICE':")

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
            self.data=pd.concat([self.data,pd.DataFrame(self.position)])
            self.data.reset_index(drop=True,inplace=True)

    def _show_position(self,verbose):
        Util.vprint(verbose,f"\tPosition Epoch: prolectic gregorian {self.epoch.datepro}, JED = {self.epoch.jed}")
        Util.vprint(verbose,f"\tCoordinates @ J2000: ")
        Util.vprint(verbose,f"\t\tEquatorial:",D2H(self.RAJ2000),D2H(self.DecJ2000))
        Util.vprint(verbose,f"\t\tEcliptic:",D2H(self.LonJ2000),D2H(self.LatJ2000))
        Util.vprint(verbose,f"\tCoordinates @ Epoch : ")
        Util.vprint(verbose,f"\t\tEquatorial:",D2H(self.RAEpoch),D2H(self.DecEpoch))
        Util.vprint(verbose,f"\t\tEcliptic:",D2H(self.LonEpoch),D2H(self.LatEpoch))
        Util.vprint(verbose,f"\tObserving conditions: ")
        Util.vprint(verbose,f"\t\tDistance to site [au]: ",self.distance)
        Util.vprint(verbose,f"\t\tDistance to sun [au]: ",self.sundistance)
        Util.vprint(verbose,f"\t\tSolar elongation [deg]: ",D2H(self.elongation))
        Util.vprint(verbose,f"\t\tPhase angle [deg]: ",D2H(self.phase))
        Util.vprint(verbose,f"\t\tMagnitude: ",self.magnitude)
        Util.vprint(verbose,f"\tOther properties: ")
        Util.vprint(verbose,f"\t\tLocal true sidereal time: ",D2H(self.site.ltst))
        Util.vprint(verbose,f"\t\tHour angle @ Epoch: ",D2H(self.HA))
        Util.vprint(verbose,f"\t\tLocal coordinates @ Epoch: ",D2H(self.az),D2H(self.el))
        
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
        self.data = pd.DataFrame()

    def _store_data(self,metstr):
        self.position.update({
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

    def _calculate_coordinates(self,tt,site):

        # Update orientation of site
        M_fixedplanet_J2000 = spy.pxform(site.planet.frameplanet,'J2000',tt)
        pos_J2000 = spy.mxv(M_fixedplanet_J2000,site.pos_fixedplanet)

        # Retrieve positions in space
        self.site_planet_SSB_J2000,lt = spy.spkezr(site.planet.id,tt,'J2000','None','SSB')
        self.planet_SSB_J2000,lt = spy.spkezr(self.id,tt,'J2000','None','SSB')
        self.sun_SSB_J2000,lt = spy.spkezr('10',tt,'J2000','None','SSB')
        self.site_SSB_J2000 = self.site_planet_SSB_J2000[:3] + pos_J2000 

        # Celestial Coordinates at J2000
        self.planet_site_J2000 = self.planet_SSB_J2000[:3] - self.site_SSB_J2000
        self.vplanet_site_J2000 = self.planet_SSB_J2000[3:] - self.site_planet_SSB_J2000[3:]
        self.planet_site_J2000 = self.planet_site_J2000
        self.vplanet_site_J2000 = self.vplanet_site_J2000

        r,RAJ2000,DECJ2000 = spy.recrad(self.planet_site_J2000)
        RAJ2000 = RAJ2000*RAD/15
        DecJ2000 = DECJ2000*RAD

        # Ecliptic coordinates J2000
        self.planet_site_EJ2000 = spy.mxv(M_J2000_ECLIPJ2000,self.planet_site_J2000)
        r,LonJ2000,LatJ2000 = spy.recrad(self.planet_site_EJ2000)
        LonJ2000 = LonJ2000*RAD
        LatJ2000 = LatJ2000*RAD

        return np.array([RAJ2000,DecJ2000,LonJ2000,LatJ2000])

    def _calculate_proper_motion(self,tt,site,dt=1*DAY,position=None):

        # Time before
        pos_t_m_dt = self._calculate_coordinates(tt-dt,site)

        # Present time
        if position is None:
            pos_t = self._calculate_coordinates(tt,site)
        else:
            pos_t = position

        # Time after
        pos_t_p_dt = self._calculate_coordinates(tt+dt,site)

        # First derivative
        dpdt = (pos_t_p_dt - pos_t_m_dt)/(2*dt)

        # Second derivative
        d2pdt2 = (pos_t_p_dt - 2*pos_t + pos_t_m_dt)/dt**2

        # Quantities are given in deg/s and deg/s^2

        return dpdt/(MARCSEC/JULYEAR),d2pdt2/(MARCSEC/JULYEAR**2)

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
        self.data = pd.DataFrame()

    def _store_data(self,metstr):
        self.position.update({
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

###############################################################
# Observing site
###############################################################
class Observer(object):
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
        self.update_site(Time(0,format='tt',scale='tt'))

    def update_site(self,mtime,verbose=False):  
        """Update site
        """

        if 'epoch' in self.__dict__.keys():
            if mtime.jed == self.epoch.jed:
                Util.vprint(verbose,f"Orientation of site not updated")
                return
            else:
                Util.vprint(verbose,f"Updating orientation of site (old time {self.epoch.datepro}, new time {mtime.datepro})")

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
        self.ltst = self.epoch.gtst + self.lon/1

###############################################################
# Sky coordinates
###############################################################
class Heka(object):


    def where_in_sky(sebau,at=None,site=None,
                     inplace=False,bar=False):
        """Compute the horizontal coordinates at an epoch of one or many celestial objects.
        
        Parameters:
            sebau: Planet | Stars | Series | DataFrame:
                Data object containing information on one or many celestial objects (sebau).

            site: montu.ObservingSite:
                Observing site.

            at: montu.Time, default = None:
                Date to which we want to precess coordinates.
                If None is because Sebau data include epoch.

            inplace: Boolean, default = False:
                If True the changes are stored in the input object (sebau).

            bar: Boolean, default = False:
                Use tqdm to show advance.

        Return:
            If inplace is False, the routine return an object of the 
            same type as `sebau` but containing the new fields:
                RAJ2000t, DecJ2000t: coordinates displaced by proper motion.
                RAEpoch, DecEpoch: coordinates precessed.
                az, el, zen: azimuth, elevation and zenith angle.

            If inplace is True, the same fields are added to sebau.
        """
        # By default we assume that the provided object is a pandas Series
        sebau_is_series = True
        # According to type of object decide the action to take
        if isinstance(sebau,Planet):
            if sebau.data is not None:
                # Recover data
                sebau = sebau.data
            else:
                # Get position among stars 
                sebau.where_among_stars(at,site)
                # Extract data 
                sebau = sebau.position
            sebau_is_series = False
        elif isinstance(sebau,Stars):
            # Extract data about the stars provided
            sebau = sebau.data
            sebau_is_series = False
            pass
        elif type(sebau) is pd.Series:
            # Check
            pass
        elif type(sebau) is pd.DataFrame:
            sebau_is_series = False
            pass
        else:
            raise AssertionError(f"Celestial object provided is not of a recognized class ({type(sebau)})")
        if sebau_is_series:
            series = sebau
            sebau = pd.DataFrame([sebau])
        if not inplace:
            sebau = copy.deepcopy(sebau)

        # Check if Sebau include epoch
        epoch = at
        epoch_implicit = False
        if 'tt' in sebau.keys():
            epoch_implicit = True
        else:
            # Update site
            site.update_site(epoch)

        # Check if sebau coordinates has been precessed
        if 'RAEpoch' not in sebau.keys():
            Heka.precess_to_epoch(sebau,at=at,inplace=True,bar=False)
        
        # Use a bar during calculations
        if bar:
            bar = tqdm.tqdm
        else:
            bar = lambda x:x

        # Main loop
        for index in bar(sebau.index):
            seba = sebau.loc[index]

            if epoch_implicit:
                # Update site if epoch implicit
                site.update_site(Time(seba.tt))

            # Get coordinates
            Dec = seba['DecEpoch']
            RA = seba['RAEpoch']

            # Compute hour angle
            HA = np.mod(site.ltst - RA,24)

            # Elevation and azimuth
            az,el = Heka._to_altaz(HA,Dec,site.lat)

            # Zenital distance
            zen = 90 - el

            # Atmospheric refraction

            # Compute airmass 
            
            # Results
            results = dict(
                HA = HA,
                az = az,
                el = el,
                zen = zen,
                epoch_tt = site.epoch.tt,
                epoch_jed = site.epoch.jed,
            )

            if not sebau_is_series:
                # Add fields to object
                sebau.loc[index,results.keys()] = results.values()
                
        # How to return results
        if not inplace:
            return sebau
        
        elif sebau_is_series:
            # Inplace for series
            for key,item in results.items():
                series[key] = item

    def precess_to_epoch(sebau,at=None,inplace=False,bar=False):
        """Precess coordinates of celestial objects (sebau in ancient egyptian)
        to epoch of observation

        Parameters:
            sebau: Series | DataFrame:
                Data containing the information on sebau.
                The object must have the following columns or attributes:
                    DecJ2000: Declination in J2000 [deg]
                    RAJ2000: Right ascension in J2000 [hours]

            mtime: montu.Time, default = None:
                Date to which we want to precess coordinates.
                If None is because Sebau data include epoch.

            inplace: Boolean, default = False:
                If True the changes are stored in the input object (sebau).

            bar: Boolean, default = False:
                Use tqdm to show advance.

        Return:
            If inplace is False, the routine return an object of the 
            same type as `sebau` but containing the new fields:
                RAJ2000t, DecJ2000t: coordinates displaced by proper motion.
                RAEpoch, DecEpoch: coordinates precessed.

            If inplace is True, the same fields are added to sebau.
        """
        # By default we assume that the provided object is a pandas Series
        sebau_is_series = True
        # According to type of object decide the action to take
        if isinstance(sebau,Planet):
            # Get position among stars 
            sebau.where_among_stars(at,site)
            # Extract data 
            sebau = sebau.position
            sebau_is_series = False
        elif isinstance(sebau,Stars):
            # Extract data about the stars provided
            sebau = sebau.data
            sebau_is_series = False
            pass
        elif type(sebau) is pd.Series:
            # Check
            pass
        elif type(sebau) is pd.DataFrame:
            sebau_is_series = False
            pass
        else:
            raise AssertionError(f"Celestial object provided is not of a recognized class ({type(sebau)})")
        if sebau_is_series:
            series = sebau
            sebau = pd.DataFrame([sebau])
        if not inplace:
            sebau = copy.deepcopy(sebau)
        
        # Check if Sebau include epoch
        epoch_implicit = False
        if 'tt' in sebau.keys():
            epoch_implicit = True

        # Check if sebau are stars
        are_stars = False
        if 'SpType' in sebau.keys():
            are_stars = True
            
        # Use a bar during calculations
        if bar:
            bar = tqdm.tqdm
        else:
            bar = lambda x:x

        # Epoch
        epoch = at

        # Main loop
        for index in bar(sebau.index):
            seba = sebau.loc[index]

            # Epoch
            if epoch_implicit:
                epoch = Time(seba.tt)

            # If stars apply motion
            if are_stars:
                dDec = seba.pmDec*(epoch.jed-JED_2000)*MARCSEC/365.25
                DecJ2000t = seba.DecJ2000 + dDec
                dRA = seba.pmRA*(epoch.jed-JED_2000)*MARCSEC/365.25/15
                RAJ2000t = seba.RAJ2000 + dRA
            else:
                # Get cordinates
                DecJ2000t = seba.DecJ2000
                RAJ2000t = seba.RAJ2000

            # Normal vector pointing to J2000 coordinates
            uJ2000 = spy.latrec(1,15*RAJ2000t*DEG,DecJ2000t*DEG)
            
            # Transform vector
            uEpoch = spy.mxv(epoch.M_J2000_Epoch,uJ2000)

            # Get spherical coordinates of new vector
            r,RA,Dec = spy.recrad(uEpoch)
            DecEpoch = Dec*RAD
            RAEpoch = RA*RAD/15

            # Get result
            results = dict(
                RAJ2000t = RAJ2000t,
                DecJ2000t = DecJ2000t,
                RAEpoch = RAEpoch,
                DecEpoch = DecEpoch,
                epoch_tt = epoch.tt,
                epoch_jed = epoch.jed,
            )

            if not sebau_is_series:
                # Add fields to object
                sebau.loc[index,results.keys()] = results.values()
                
        if not inplace:
            return sebau
        
        elif sebau_is_series:
            for key,item in results.items():
                series[key] = item
    
    def _to_altaz(HA,Dec,lat):
        """Transform from local equatorial to azimuth and elevation
        """
        # Elevation and azimuth
        el = np.arcsin(np.sin(Dec*DEG)*np.sin(lat*DEG) + \
                    np.cos(Dec*DEG)*np.cos(lat*DEG)*np.cos(HA*15*DEG))*RAD
        az = np.arctan2(-np.sin(HA*15*DEG)*np.cos(Dec*DEG)/np.cos(el*DEG),
                        (np.sin(Dec*DEG) - np.sin(lat*DEG)*np.sin(el*DEG))/\
                            (np.cos(lat*DEG)*np.cos(el*DEG)))*RAD
        az = np.mod(az,360)
        return az,el
    
    def _daily_motion(HA0,Dec0,lat,deltat,tai2sid=TAI_TO_SID,
                      pm_ra=0,pm_dec=0,pa_ra=0,pa_dec=0):
        """Propagate local position of a body at a given sky coordinates

        Parameters:
            HA0 [hours], Dec0 [deg]:
                Initial local sky coordinates of the object.

            lat: float [deg]:
                Latitude of the observing site

            pm_ra, pm_dec: float [uarcsec], default = 0:
                Propet motion in RA and in Dec.

            tai2sid: float, default = TAI_TO_SID:
                Ratio of sidereal day to TAI day.
                It can be computed using Time.tai_to_sid(mtime).

        Return:
            HA: float [hour]:
                Hour angle.
            az, el: float [degree]:
                Azimuth and elevation            
        """
        # Uodate proper motion (leap frog)
        pm_ra = pm_ra + pa_ra*(deltat/JULYEAR)
        pm_dec = pm_dec + pa_dec*(deltat/JULYEAR)

        # Equatorial coordinates
        Dec0 = Dec0 + deltat*(pm_dec*MARCSEC/JULYEAR)
        HA0 = HA0 + deltat*(pm_ra*MARCSEC/JULYEAR)/15

        # Advance in sidereal time
        dst = tai2sid*deltat/HOUR # Advanced hours
        
        # Hour angle after deltat
        HA = np.mod(HA0 + dst,24)
        
        # Elevation and azimuth
        az,el = Heka._to_altaz(HA,Dec0,lat)
        
        return HA,az,el,pm_ra,pm_dec

    def move_over_nut(seba,at=None,site=None,
                      during=24*HOUR,each=1*HOUR,
                      tai2sid=None):
        """Propagate local position of a body at a given sky coordinates.

        Parameters:
            seba: pandas.Series or pandas.DataFrame.
                Object containing the coordinates of the celestial object.
                The object should have as attributes: RAEpoch, DecEpoch, pmRA and pmDec.

            site: montu.ObservingSite:
                Observing site.

            during: float [seconds], default = 24*HOUR:
                during of the motion in seconds.

            each: float [seconds], default = 1*HOUR
                Step size of the motion.

            tai2sid: float, default = None:
                Ratio of sidereal day to TAI day.
                If None the routine automatically compute it using 
                Time.tai_to_sid(mtime).
        
        Return:

            path: pandas.DataFrame:
                A DataFrame containing the following columns:
                    tt: Terrestrial time.
                    jed: Julian day in UTC scale.
                    HA: Hour angle at tt.
                    az, el, zen: Azimuth and elevation at tt
                    isvisible: True if object is reachable (is above horizon).
        """

        # Check if seba is DataFrame:
        if type(seba) == pd.DataFrame:
            if len(seba)>1:
                raise AssertionError(f"Object provided has more than 1 position ({len(seba)})")    
            # Convert seba to pandas Series
            seba = seba.iloc[0]

        elif type(seba) == pd.Series:
            # Default format
            pass

        else:
            raise AssertionError("Object provided is not DataFrame or Series")

        # Update site
        site.update_site(at)
        if tai2sid is None:
            tai2sid = Time.tai_to_sid(at)

        # Initial time
        tt0 = at.tt

        # Check if sky coordinates has been precessed
        if 'RAEpoch' in seba.keys():
            dt = abs(tt0 - seba.epoch_tt)
            if dt>1:
                raise ValueError(f"The epoch provided JD {at.jed} does not coincide with epoch of sky coordinates {seba.epoch_jed}")
        else:
            Heka.precess_to_epoch(seba,at=at,inplace=True)

        # Get coordinates
        Dec0 = float(seba.DecEpoch)
        RA0 = float(seba.RAEpoch)
        HA0 = np.mod(site.ltst - RA0,24)

        # Get proper motion and acceleration
        if 'pmRA' in seba.keys():
            pmRA = seba.pmRA
            pmDec = seba.pmDec
        if 'paRA' in seba.keys():
            paRA = seba.paRA
            paDec = seba.paDec

        # Main loop
        results = []
        for i,tt in enumerate(Util.arange(tt0,tt0+during,each)):

            HA,az,el,pmRA,pmDec = Heka._daily_motion(HA0,Dec0,site.lat,tt-tt0,
                                                     tai2sid=tai2sid,pm_ra=pmRA,pm_dec=pmDec)
            zen = 90 - el

            # Results
            results += [dict(
                tt = tt,
                HA = HA,
                az = az,
                el = el,
                zen = zen,
                isvisible = True if el>0 else False,
            )]
            
        path = pd.DataFrame(results)
        return path

    def plot_over_nut(paths):
        
        if not isinstance(paths,list):
            paths = [paths]
        
        # Create plot
        fig,axs = plt.subplots(subplot_kw=dict(projection='polar'),
                                facecolor='g')
        axs.set_facecolor('black')

        # Plot paths
        for path in paths:
            # Choose only points above horizon (reachable)
            cond = (path.el >= 0)
            # Plot 
            axs.scatter(path[cond].az*DEG,path[cond].zen,
                        color='y') #,s=path[cond].el)

        # Decoration
        # Change from zenithal angle to elevation labels
        axs.set_rgrids(np.arange(0,90,15),angle=0)
        el_labels = []
        zs = axs.get_yticks()
        for z in zs:
            el_labels += [f'+{90-z}']
        axs.set_yticklabels(el_labels,color='w',fontsize=6);
        
        # Change azimut labels
        azs = np.unique(list(np.arange(0,360,15))+[0,90,180,270])*DEG
        axs.set_xticks(azs)
        azlabels = []
        for az in azs:
            az = az*RAD
            if az==0:azlabel='N'
            elif az==90:azlabel='E'
            elif az==180:azlabel='S'
            elif az==270:azlabel='W'
            else:azlabel=f'{az:.0f}'
            azlabels += [azlabel]
        axs.set_xticklabels(azlabels,fontsize=6,color='w');

        # Grid
        axs.grid(alpha=0.3)
        
        return fig,axs