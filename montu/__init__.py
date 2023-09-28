import numpy as np
import spiceypy as spy
import os
import requests
from tqdm import tqdm 
import pandas as pd
import matplotlib.pyplot as plt
import copy
from astropy.time import Time

# Constants
RAD = 180/np.pi
DEG = 1/RAD

DAY = 86400

# Version
from montu.version import *

# Global variables
KERNELS = {
    'latest_leapseconds.tls':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls',
    'de441_part-1.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp',
    'de441_part-2.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp',
}
KERNELS_UP = dict()

class Montu(object):
    # Load kernels
    def load_kernels(dir='/tmp/',kernels=KERNELS,verbose=True):
        # Check if dir exists
        if not os.path.exists(dir):
            os.system(f"mkdir -p {dir}")    
        for kernel,item in kernels.items():
            kernel_path = dir+"/"+kernel
            if not os.path.exists(kernel_path):
                if verbose:print(f"Downloading '{kernel}'...")
                Montu._wget(item,kernel_path)
            else:
                if kernel not in KERNELS_UP.keys():
                    if verbose:print(f"Loading kernel {kernel}")
                    spy.furnsh(kernel_path)
                    KERNELS_UP[kernel] = True

    # Util routines
    def _wget(url, file_name):
        """
        Source: ChatGPT
        """
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))

        # Initialize the progress bar
        progress_bar = tqdm(total=total_size, unit='B', unit_scale=True)

        with open(file_name, 'wb') as file:
            for data in response.iter_content(chunk_size=1024):
                file.write(data)
                progress_bar.update(len(data))

        progress_bar.close()

    def _data_path(filename):
        """
        Get the full path of the `datafile` which is one of the datafiles provided with the package.
        
        Parameters:
            filename: Name of the data file, string.
            
        Return:
            Full path to package datafile in the python environment.
            
        """
        return os.path.join(os.path.dirname(__file__),'data',filename);

    def dec2hex(dec):
        sgn = np.sign(dec)
        dec = abs(dec)
        h = int(dec)
        mf = 60*(dec - int(dec))
        m = int(mf)
        s = 60*(mf - m)
        return sgn*h,m,s

    def haversine_distance(lat1, lon1, lat2, lon2):
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        return c
    
    def dt2cal(dt,bce=False):
        """Convert array of datetime64 to a calendar array of year, month, day, hour,
        minute, seconds, microsecond with these quantites indexed on the last axis.

        Parameters
        ----------
        dt : datetime64 array (...)
            numpy.ndarray of datetimes of arbitrary shape

        Returns
        -------
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

        out = np.array([float(o) for o in out])
        out[0] = -out[0] if bce else out[0]
        return out

class Stars(object):

    def __init__(self,data=None,filename=None):
        if data is not None:
            self.data = copy.deepcopy(data)
        elif filename:
            self.data = pd.read_csv(filename)
        else:
            self.data = pd.read_csv(Montu._data_path('stars.csv'))

    def get_stars(self,**args):
        if len(args)==0:
            return self.data
        cond = np.array([True]*len(self.data))
        for key,item in args.items():
            # print("Antes: ",cond.sum())
            if isinstance(item,list):
                min = float(item[0])
                max = float(item[1])
                cond = (self.data[key]>=min)&(self.data[key]<=max)&(cond)
            else:   
                cond = (self.data[key]==item)&(cond)
            # print("Después: ",cond.sum())
        return Stars(self.data[cond])
    
    def get_stars_area(self,RA=0,Dec=0,radius=10,**kwargs):
        stars = self.get_stars(
            RA=[RA-radius/15,RA+radius/15],
            Dec=[Dec-radius,Dec+radius],
            **kwargs
        )
        return stars
    
    def plot_stars(self,pad=0):
        fig,axs = plt.subplots(1,1,figsize=(5,5))
        axs.set_facecolor('black')
        axs.scatter(15*self.data.RA,self.data.Dec,marker='*',color='y',s=10)

        # Decoration
        axs.grid(alpha=0.2)
        axs.tick_params(axis='x',direction='in', pad=-10, size=10)
        axs.tick_params(axis='y',direction='in', pad=-20, size=10)
        for axis in axs.xaxis,axs.yaxis:
            for t in axis.get_ticklabels():
                t.set_color('w')
                t.set_fontsize(6)

        axs.set_xlabel('RA [deg]',labelpad=-20,color='w',fontsize=6)
        axs.set_ylabel('Dec [deg]',labelpad=-30,color='w',fontsize=6)
        axs.axis('equal')

        rang = max(15*((self.data.RA).max()-(self.data.RA).min()),
                (self.data.Dec).max()-(self.data.Dec).min())
        axs.margins(pad*rang)
        fig.tight_layout()
        return fig,axs

class MonTime(object):
    
    def __init__(self,datestr,format='iso'):
        self._parse_datestr(datestr,format)
    
    def _parse_datestr(self,datestr,format):
        """Create a time object with a date in string

        Parameters
            datestr: string:
                String in format [-]CCYYY-MM-DD HH:MM:SS.sss
        
        Updated attributes:
            datestr: 
            bce: is date bce (before current era or befor Christ)
            datetime64: numpy datetime
            cal: calendar items (year, month, day, hour, minutes, seconds, miliseconds)
            datespice: date in SPICE format
            et: ephemeris time [seconds]
            deltet: delta-time for the epoch.
            jd: julian day (TDB)
        """
        # Strip blank spaces
        self.datestr = datestr.strip()

        # Strip negative values in the year
        self.bce = True if self.datestr[0] == '-' else False
    
        # Convert to datetime64
        self.datetime64 = np.datetime64(self.datestr)

        # Extract calendar components
        self.cal = Montu.dt2cal(np.datetime64(self.datestr.strip('-')),bce=self.bce)
        
        # Write 
        suffix = f"{int(self.cal[1]):02d}-{int(self.cal[2]):02d} {int(self.cal[3]):02d}:{int(self.cal[4]):02d}:{int(self.cal[5]):02d}.{int(self.cal[6]):d}"
        if self.bce:
            self.datespice = f"{int(abs(self.cal[0])+1):d} B.C. {suffix}"
        else:
            self.datespice = f"{int(abs(self.cal[0])):d}-{suffix}"

        # Load leap-seconds kernel if not up
        if 'latest_leapseconds.tls' not in KERNELS_UP.keys():
            kernels = {key: KERNELS[key] for key in KERNELS.keys() & {'latest_leapseconds.tls'}}
            Montu.load_kernels(kernels=kernels,verbose=False)

        # Convert to et
        self.et = spy.utc2et(self.datespice)
        self.deltat = spy.deltet(self.et,"ET")
        self.et -= self.deltat
        self.jd = spy.unitim(self.et,"ET","JDTDB")

        # Astropy
        self.astrotime = Time(self.jd,format='jd',scale='tdb')
            
class Planet(object):
    def __init__(self,id):
        self.id = id

    def calc_ephemerides(self,epochs=0,location='399'):
        # Calc state vector
        X,lt = spy.spkezr('4',epochs,'J2000','LT+S',location)
        
        # Compute RA and DEC (J2000)
        r,self.RAJ2000,self.DECJ2000 = spy.reclat(X[:3])
        self.RAJ2000 = np.mod(self.RAJ2000*RAD,360)/15
        self.DECJ2000 *= RAD 

        return
        # Compute ecliptic coordinates (J2000)
        Mj2000_to_eclipJ2000 = spy.pxform('J2000','ECLIPJ2000',epochs) 
        poseclip = spy.mxv(Mj2000_to_eclipJ2000,X[:3])
        r,self.lambJ2000,self.betaJ2000 = spy.reclat(poseclip)
        self.lambJ2000 = np.mod(self.lambJ2000*RAD,360)
        self.betaJ2000 *= RAD

        # Compute horizontal coordinates for the epochs
