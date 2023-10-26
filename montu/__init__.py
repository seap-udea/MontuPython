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
# Import stars modules
###############################################################
from montu.util import *
from montu.time import *
from montu.stars import *
from montu.planets import *
from montu.observer import *
from montu.heka import *

###############################################################
# Data initialization
###############################################################
# Correction for JED
jed_correction_data = np.loadtxt(Util._data_path('corrections_dt.dat'))
JED_CORRECTION = interp1d(jed_correction_data[:,0],jed_correction_data[:,1])

# Load basic kernels
Util.load_kernels()

###############################################################
# Tests 
###############################################################
if __name__ == '__main__':
    print('Montu Python Test Suite')
