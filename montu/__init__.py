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
from pyplanets.core.coordinates import mean_obliquity as pyplanets_obliquity
from astropy import constants as astropy_constants
from astropy import units as monu # Montu Units are astropy units
from astropy.time import Time as astropy_Time
from astroquery.jplhorizons import Horizons as astropy_Horizons

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
            
            else:
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
    
    def dec2hex(dec):
        sgn = np.sign(dec)
        dec = abs(dec)
        h = int(dec)
        mf = 60*(dec - int(dec))
        m = int(mf)
        s = 60*(mf - m)
        return sgn*h,m,s

###############################################################
# MonTime Class
###############################################################
class MonTime(object):
    """Create a time-frame object
    
    This is one of the most important classes in the package, since
    it manages times and dates.

    Initialization parameters:
        date: string | float
            Date

        format: string, default = 'iso'
            Format of the input date. Possible values: 'iso', 'tt', 'jd'
        
        scale: string, default = 'utc'
            Scale of the time: 'tt' (terrestrial time, uniform), 
            'utc' (coordinated, based on rotation)

        proleptic: boolean, default = True
            Proleptic gregorian correspond to the case when the Gregorian calendar
            is extended before the adoption date at 1582-10-15. 

            When proleptic = False, the Julian calendar is used before the adoption
            date.

    Attributes:
        et: float
            Ephemeris time in tt. 
        deltat: float
            Difference Dt = TT - UTC. 
        jdt: float
            Julian day in terrestrial time.
        jde: float
            Julian day in UTC.
        dateastro: string
            Date in astropy format '[-]CCYY-MM-DD HH:MM:SS.fff'
            
        datestr: string
            Date in gregorian prolectic, with format '[bce]CCYY-MM-DD HH:MM:SS.fff'
        datemix: string
            Date in a mixed style (non-prolectic), with format '[bce]CCYY-MM-DD HH:MM:SS.fff', 
            meaning that when date is previous to 1582-10-15 the date is given in Julian 
            Calendar.
        dateastro: string
            Date in gregorian prolectic but in astronomical format '[-]CCYY-MM-DD HH:MM:SS.fff'

    Other attributes (related to frame of reference)
        obliquity: float
            Mean obliquity for the date.
        M_J2000_Epoch: float
            Matrix transforming from J2000 to Epoch.
            
    Examples:

        Initialization using a string:
        
            MonTime('-2500-01-01 12:00:00.00',format='iso',scale='utc',proleptic=True)
            MonTime('bce 2501-01-01 12:00:00.00',format='iso',scale='utc',proleptic=True)
            MonTime('2501 B.C. 01-01 12:00:00.00',format='iso',scale='tt',proleptic=True)
            MonTime('-2500-01-01 12:00:00.00',format='iso',scale='utc',proleptic=True)
            MonTime('-2500-01-22 12:00:00.00',format='iso',scale='utc',proleptic=False)

        Initialization using a float:

            MonTime(0,format='tt',scale='tt') # Gives 2000-01-01 12:00:00 TT 
            MonTime(2451545,format='jd',scale='tt') # Gives 2000-01-01 12:00:00 TT 
    """    
    def _parse_datestr(self,date,verbose=False):
        """Parse string
        """
        
        # Default format
        self.style = 'astro' # Default style
        self.dateastro = date # Default date

        # Strip blank spaces
        self.date = date.strip()

        # Is time before current era
        self.bce = False
        if self.date[0] == '-':
            self.bce = True
        elif 'b' in self.date.lower():
            self.bce = True
            self.style = 'calendar'

        # Convert all formats to dateastro
        if self.bce and (self.style == 'calendar'):
            subs1 = lambda m:str(-(int(m.group(1))-1))
            subs2 = lambda m:str(-(int(m.group(1))-1))+'-'
            self.dateastro = re.sub('bc[a-z\s]*(\d+)',subs1,self.dateastro.lower())
            self.dateastro = re.sub('(\d+)\s*b[\.]*c[\.]*\s*',subs2,self.dateastro.lower())

        Montu.vprint(verbose,"Date astro:",self.dateastro)

        # Create calendar and datetime object
        self.datetime64 = np.datetime64(self.dateastro)
        Montu.vprint(verbose,"Datetime64:",self.datetime64)

        self.calendar = Montu.dt2cal(np.datetime64(self.dateastro.strip('-')),bce=self.bce)
        Montu.vprint(verbose,"Calendar:",self.calendar)

        self.datetime = datetime(*self.calendar[1:])
        Montu.vprint(verbose,"Datetime:",self.datetime)

        # SPICE format
        if self.bce:
            datelist = [int(f) for f in self.datetime.strftime('%Y,%m,%d,%H,%M,%S,%f').split(',')]
            datelist[0] += 1
            dateprev = datetime(*tuple(datelist))
        else:
            dateprev = self.datetime
        self.datespice = dateprev.strftime('%Y B.C. %m-%d %H:%M:%S.%f') if self.bce else self.dateastro
        Montu.vprint(verbose,"Date SPICE:",self.datespice)
        
    def __init__(self,date,format='iso',scale='tt',proleptic=True,verbose=False):
        
        self.scale = scale
        self.verbose = verbose

        if format == 'iso':

            if proleptic:
                # If proleptic use SPICE
                Montu.vprint(self.verbose,f"Proleptic {self.scale}")
                
                # Parse string
                self._parse_datestr(date,verbose=self.verbose)

                # Calculated deltat
                self.deltat = pyplanets_Epoch.tt2ut(self.calendar[0]*self.calendar[1],self.calendar[2])

                # Convert to TT
                et = spy.str2et(self.datespice)
                et -= spy.deltet(et,'ET')
                
                # Get JED
                if scale == 'tt':
                    self.tt = et
                    self.et = et - self.deltat
                else:
                    self.tt = et + self.deltat
                    self.et = et
                self.jed = spy.unitim(self.et,'ET','JED')
                self.jtd = spy.unitim(self.tt,'ET','JED')

                # Get mixed epoch
                self.mixed = pyplanets_Epoch(self.jed)
                
            else:
                # If not proleptic (mixed gregorian-julian) use PyPlanets
                Montu.vprint(self.verbose,f"Non-proleptic {self.scale}")

                # Parse string
                self._parse_datestr(date,verbose=False)

                # Calculated deltat
                self.deltat = pyplanets_Epoch.tt2ut(self.calendar[0]*self.calendar[1],self.calendar[2])

                # Create a mixed epoch
                args = (self.calendar[0]*self.calendar[1],
                        self.calendar[2],
                        self.calendar[3]+\
                        self.calendar[4]/24+\
                        self.calendar[5]/(24*60)+\
                        self.calendar[6]/(24*60*60)+\
                        self.calendar[7]/(24*60*60*1e6))
                self.mixed = pyplanets_Epoch(*args)
                
                # According to scale
                jd = self.mixed.jde()
                et = spy.unitim(jd,'JED','ET')
                
                # Choose according to scale
                if scale == 'tt':
                    self.tt = et
                    self.jtd = jd
                    self.et = self.tt - self.deltat
                    self.jed = self.jtd - self.deltat/DAY
                else:
                    self.et = et
                    self.jed = jd
                    self.tt = self.et + self.deltat
                    self.jtd = self.jed + self.deltat/DAY
                
                c = spy.et2utc(et+spy.deltet(et,'ET'),'C',4)
                            
                # Use the original et without modification
                if self.bce:
                    d = datetime.strptime(c,'%Y B.C. %b %d %H:%M:%S.%f')
                    datestr = f'{-(d.year-1)}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}.{d.microsecond}'
                    
                else:
                    d = datetime.strptime(c,'%Y %b %d %H:%M:%S.%f')
                    datestr = f'{d.year}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}.{d.microsecond}'
    
                self._parse_datestr(datestr,verbose=self.verbose)

            Montu.vprint(self.verbose,"Mixed:",self.mixed.get_full_date())
            Montu.vprint(self.verbose,"et = ",self.et)
            Montu.vprint(self.verbose,"JED = ",self.jed)
            Montu.vprint(self.verbose,"Delta-t:",self.deltat)
            Montu.vprint(self.verbose,"tt = ",self.tt)
            Montu.vprint(self.verbose,"JTD = ",self.jtd)

        else:
            # Initialize dates using terrestrial time
            if format == 'jd':
                jd = date
            elif format == 'tt':
                et = date
                jd = spy.unitim(et,'ET','JED')
            else:
                raise AssertionError(f"Format '{format}' not recognized (valid 'iso', 'tt', 'jd')")

            # Initialize Epoch
            self.mixed = pyplanets_Epoch(jd)
            year,month,day = self.mixed.get_date()
            self.bce = True if year<=0 else False

            self.deltat = pyplanets_Epoch.tt2ut(year,month)

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

            # Date string
            c = spy.et2utc(et+spy.deltet(et,'ET'),'C',4)
                        
            # Use the original et without modification
            if self.bce:
                d = datetime.strptime(c,'%Y B.C. %b %d %H:%M:%S.%f')
                datestr = f'{-(d.year-1)}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}.{d.microsecond}'
        
            else:
                d = datetime.strptime(c,'%Y %b %d %H:%M:%S.%f')
                datestr = f'{d.year}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}.{d.microsecond}'

            self._parse_datestr(datestr,verbose=self.verbose)

            Montu.vprint(self.verbose,"Mixed:",self.mixed.get_full_date())
            Montu.vprint(self.verbose,"et = ",self.et)
            Montu.vprint(self.verbose,"JED = ",self.jed)
            Montu.vprint(self.verbose,"Delta-t:",self.deltat)
            Montu.vprint(self.verbose,"tt = ",self.tt)
            Montu.vprint(self.verbose,"JTD = ",self.jtd)

        # Create astrotime 
        self.astrotime = astropy_Time(self.jtd,format='jd',scale='tdb')

        # Obliquity
        self.epsilon = pyplanets_obliquity(self.mixed)

        # Update matrices
        self.M_J2000_Epoch = spy.pxform('J2000','EARTHTRUEEPOCH',self.tt)
        
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
        axs.set_xlabel('RA [deg]',fontsize=10)
        axs.set_ylabel('Dec [deg]',fontsize=10)
        
        # Range
        rang = max(15*((self.data.RA).max()-(self.data.RA).min()),
                (self.data.Dec).max()-(self.data.Dec).min())
        axs.margins(pad*rang)
        
        axs.grid(alpha=0.2)
        axs.axis('equal')
        fig.tight_layout()
        return fig,axs


###############################################################
# Planetary body
###############################################################
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