import numpy as np
import spiceypy as spy
import os
import requests
from tqdm import tqdm 
import pandas as pd
import matplotlib.pyplot as plt
import copy
from astropy import constants
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, FK5, ICRS
import astropy.units as u
from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator
from functools import lru_cache
import warnings
warnings.filterwarnings("ignore")

LEGACY = True
MAIN = False

# Constants
RAD = 180/np.pi
DEG = 1/RAD

DAY = 86400 # s
YEAR = 365.25*DAY # s
KM = 1e3 # m

AU = 149597870.700 # km
CSPEED = 299792.458 # km/s 

# Version
from montu.version import *

# Global variables
KERNELS = {
    'latest_leapseconds.tls':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls',
    'de441_part-1.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp',
    'de441_part-2.bsp':'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp',
    'frame.tk':'',
    'pck00011.tpc':'',
    'earth_assoc_itrf93.tf':''
}
KERNELS_UP = dict()

class Montu(object):

    # Print with verbosity
    def vprint(verbose=True,*args):
        if verbose:print(*args)

    # Load kernels
    def load_kernels(dir='/tmp/',kernels=KERNELS,verbose=True):
        # Check if dir exists
        if not os.path.exists(dir):
            os.system(f"mkdir -p {dir}")    

        # Load kernel
        for kernel,item in kernels.items():
            if len(item) == 0:
                if verbose:print(f"Loading local kernel {kernel}")
                spy.furnsh(Montu._data_path(kernel))
                continue
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
    
    def time_ticks(ax,nticks=None,**kwargs):
        if nticks:
            etmin, etmax = ax.get_xlim()
            eticks = np.linspace(etmin,etmax,nticks)
            ax.set_xticks(eticks)
        else:
            eticks = ax.get_xticks()
        elabels = []
        for etick in eticks:
            mt = MonTime(etick,format='et')
            elabels += [f"{int(abs(mt.cal[0])):d}-{int(mt.cal[1]):02d}-{int(mt.cal[2]):02d}"]
        ax.set_xticklabels(elabels,**kwargs)

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
    
    def plot_stars(self,pad=0,**kwargs):

        # Default options
        default = dict(figsize=(5,5))
        default.update(kwargs)

        fig,axs = plt.subplots(1,1,**default)
        axs.set_facecolor('black')
        axs.scatter(15*self.data.RA,self.data.Dec,marker='*',color='y',s=10)

        # Decoration
        axs.grid(alpha=0.2)
        axs.tick_params(axis='x',direction='in', pad=-10, size=10)
        axs.tick_params(axis='y',direction='in', pad=-20, size=10)
        for axis in axs.xaxis,axs.yaxis:
            for t in axis.get_ticklabels():
                t.set_color('w')
                t.set_fontsize(10)

        axs.set_xlabel('RA [deg]',labelpad=-20,color='w',fontsize=10)
        axs.set_ylabel('Dec [deg]',labelpad=-30,color='w',fontsize=10)
        axs.axis('equal')

        rang = max(15*((self.data.RA).max()-(self.data.RA).min()),
                (self.data.Dec).max()-(self.data.Dec).min())
        axs.margins(pad*rang)
        fig.tight_layout()
        return fig,axs

class MonTime(object):
    
    def __init__(self,date,format='iso',scale='utc'):
        if format == 'iso':
            self._parse_datestr(date,format,scale)
        elif format == 'spice':
            self._parse_datestr(date,format,scale)
        elif format == 'et':
            self._convert_from_unitim(date,format,scale)
        elif format == 'jd':
            self._convert_from_unitim(date,format,scale)
        else:
            raise AssertionError(f"Format '{format}' not supported")
    
    def _convert_from_unitim(self,utime,format='et',scale='utc'):
        if format == 'jd':
            # Julian day
            utime = spy.unitim(utime,"JDTDB","ET")

        # Ephemeris time
        self.et = utime
        self.deltat = spy.deltet(utime,"ET")
        if scale == 'utc':
            self.et -= self.deltat
        
        # SPICE format
        self.datespice = spy.et2utc(self.et+self.deltat,'C',3)

        # Julian day
        self.jd = spy.unitim(self.et,"ET","JDTDB")

        # Astropy
        self.astrotime = Time(self.jd,format='jd',scale='tdb')
        self.datestr = self.astrotime.iso

        # Is time before current era
        self.bce = True if self.datestr[0] == '-' else False

        # Convert to datetime64
        self.datetime64 = np.datetime64(self.datestr)

        # Extract calendar components
        self.cal = Montu.dt2cal(np.datetime64(self.datestr.strip('-')),bce=self.bce)
    
    def _parse_datestr(self,datestr,format='iso',scale='utc'):
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

        # Is time before current era
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
        if scale == 'tt':
            self.et -= self.deltat
        
        self.jd = spy.unitim(self.et,"ET","JDTDB")

        # Astropy
        self.astrotime = Time(self.jd,format='jd',scale='tdb')
            
class Planet(object):
    def __init__(self,id):
        self.id = id

    def calc_ephemerides(self,epochs=0,location='399'):
        # Calc state vector
        self.X,lt = spy.spkezr(self.id, epochs,'J2000','LT+S',location)
        
        # Compute RA and DEC (J2000)
        r,self.RAJ2000,self.DECJ2000 = spy.reclat(self.X[:3])
        self.RAJ2000 = np.mod(self.RAJ2000*RAD,360)/15
        self.DECJ2000 *= RAD 

        # Compute ecliptic coordinates (J2000)
        Mj2000_to_eclipJ2000 = spy.pxform('J2000','ECLIPJ2000',epochs) 
        poseclip = spy.mxv(Mj2000_to_eclipJ2000,self.X[:3])
        r,self.lambJ2000,self.betaJ2000 = spy.reclat(poseclip)
        self.lambJ2000 = np.mod(self.lambJ2000*RAD,360)
        self.betaJ2000 *= RAD

        # Compute horizontal coordinates for the epochs

    def angle_to_planet(self,planet):
        # Get position vector
        rself = self.X[:3]
        rplan = planet.X[:3]

        # Get cosine of angle
        costeta = (rself@rplan)/((rself@rself)**0.5*(rplan@rplan)**0.5)
        return np.arccos(costeta)*RAD
    
    def angle_to_star(self,star):
        # Calculate the angular distance
        angdist = Montu.haversine_distance(self.DECJ2000*DEG,15*self.RAJ2000*DEG,
                                           float(star.data.Dec)*DEG,15*float(star.data.RA)*DEG)*RAD
        return angdist

    def calc_proper_motion(self,et,dt=1):

        # Get consecutive positions
        self.calc_ephemerides(et-dt)
        ra_ini = 15*self.RAJ2000
        dec_ini = self.DECJ2000
        self.calc_ephemerides(et+dt)
        ra_end = 15*self.RAJ2000
        dec_end = self.DECJ2000

        # Compute distance traversed 
        angdist = Montu.haversine_distance(dec_ini*DEG,ra_ini*DEG,dec_end*DEG,ra_end*DEG)*RAD
        
        # Compute angular velocities
        ra_speed = (ra_end - ra_ini)/(2*dt)
        dec_speed = (dec_end - dec_ini)/(2*dt)
        ang_speed = angdist/(2*dt)

        return ra_speed,dec_speed,ang_speed

class Consts(object):
    pass

Consts.au = constants.au.value

class PlanetaryBody(object):
    """Create a planetary body

    Examples: 
        earth = PlanetaryBody(id='399')
        print(Earth.f)
    """

    @lru_cache()
    def __new__(cls,id):
        """This method is intended to avoid creating a new object with the same id
        Instead this method create a clone of the previously created object.
        """
        return super().__new__(cls)

    def __init__(self,id,framebody='IAU_EARTH'):

        # Obtain properties of object
        self.id = id
        n,rs=spy.bodvrd(self.id,"RADII",3)
        self.Re=rs[0]
        self.Rt=rs[1]
        self.Rp=rs[2]
        self.f=(self.Re-self.Rp)/self.Re
        self.framebody = framebody

        # Set dummy ephemerides
        self.epochs = []
        self.eq_J2000 = SkyCoord(ra=0*u.deg,dec=0*u.deg)
        self.eq_epoch = SkyCoord(ra=0*u.deg,dec=0*u.deg)
        self.az_epoch = None

    def calculate_ephemerides(self,location,epochs,verbose=True):

        self.epochs = epochs

        # Ephemeris time
        et = epochs.et 
        
        # Retrieve position of planet
        try:
            Montu.vprint(verbose,f"Retrieving position for object {location.planet.id}")
            planetSSBJ2000,lt = spy.spkezr(location.planet.id,et,'J2000','None','SSB')
        except:
            Montu.vprint(verbose,f"\tCorrecting id to {location.planet.id[0]}")
            planetSSBJ2000,lt = spy.spkezr(location.planet.id[0],et,'J2000','None','SSB')
            
        try:
            Montu.vprint(verbose,f"Retrieving position for object {self.id}")
            bodySSBJ2000,lt = spy.spkezr(self.id,et,'J2000','None','SSB')
        except:
            Montu.vprint(verbose,f"\tCorrecting id to {self.id[0]}")
            bodySSBJ2000,lt = spy.spkezr(self.id[0],et,'J2000','None','SSB')

        

        """
        This has been tested using Horizons @ 2000-01-01 12:00:00.0 UTC and it gives the 
        position exactly 
        """
        M_J2000_ECLIPJ2000 = spy.pxform('J2000','ECLIPJ2000',et)
        
        Montu.vprint(verbose,f"Planet {location.planet.id} position w.r.t. SSB = ",
                     spy.mxv(M_J2000_ECLIPJ2000,planetSSBJ2000[:3]),
                     spy.mxv(M_J2000_ECLIPJ2000,planetSSBJ2000[3:])
                     )
        Montu.vprint(verbose,f"Body {self.id} position w.r.t. SSB = ",
                     spy.mxv(M_J2000_ECLIPJ2000,bodySSBJ2000[:3]),
                     spy.mxv(M_J2000_ECLIPJ2000,bodySSBJ2000[3:]),
                     )
        self.P_SSBJ2000 = planetSSBJ2000
        self.X_SSBJ2000 = bodySSBJ2000

        # Position of the observing site with respect to SSB
        location.SSBJ2000 = planetSSBJ2000[:3] + location.J2000
        
        # Celestial Coordinates at J2000
        bodyTOPOJ2000 = bodySSBJ2000[:3] - location.SSBJ2000
        r,RAbodyJ2000,DECbodyJ2000 = spy.recrad(bodyTOPOJ2000)
        self.eq_J2000 = SkyCoord(ra=RAbodyJ2000*u.rad,dec=DECbodyJ2000*u.rad)
        Montu.vprint(verbose,"Coordinates @ J2000: ",self.eq_J2000.ra/15,self.eq_J2000.dec)

        # Celestial Coordinates at Epoch
        if location.planet.id == '399':
            # Method 1: using astropy
            self.eq_epoch = self.eq_J2000.transform_to(FK5(equinox=self.epochs.astrotime))
            #"""
            # There is a second method to precess the coordinates using SPICE
            Montu.vprint(verbose,f"Coordinates @ Epoch using astropy : ",self.eq_epoch.ra/15,self.eq_epoch.dec)
            # Method 2: using SPICE
            M_J2000_Epoch = spy.pxform('J2000','EARTHTRUEEPOCH',et)
            bodyTOPOEpoch = spy.mxv(M_J2000_Epoch,bodyTOPOJ2000)
            r,RAbody,DECbody = spy.recrad(bodyTOPOEpoch)
            self.eq_epoch = SkyCoord(ra=RAbody*u.rad,dec=DECbody*u.rad)
            Montu.vprint(verbose,f"Coordinates @ Epoch using SPICE : ",self.eq_epoch.ra/15,self.eq_epoch.dec)
            #"""
        else:
            self.eq_epoch = SkyCoord(ra=RAbodyJ2000*u.rad,dec=DECbodyJ2000*u.rad)
        
        Montu.vprint(verbose,"Coordinates @ Epoch: ",self.eq_epoch.ra/15,self.eq_epoch.dec)

        # AltAz Coordinates at body
        self.az_epoch = self.eq_J2000.transform_to(location.frame)
        Montu.vprint(verbose,"AltAz @ Epoch: ",self.az_epoch.az.value,self.az_epoch.alt.value)
        Montu.vprint(verbose,"AltAz @ Epoch: ",Montu.dec2hex(self.az_epoch.az.value),Montu.dec2hex(self.az_epoch.alt.value))

        # AltAz Coordinates at body (using SPICE)
        location.update_site(self.epochs)
        az,el = location.to_local(self.eq_epoch.ra.value,self.eq_epoch.dec.value)
        Montu.vprint(verbose,"AltAz @ Epoch (SPICE): ",Montu.dec2hex(az),Montu.dec2hex(el))
    
    def __str__(self):
        str = f"""
Planetary body id: {self.id}
Size: {self.Re,self.Rt,self.Rp,self.f}
Epoch: {self.epochs.datestr}
Coordinates:
    equatorial J2000 = {self.eq_J2000}
    equatorial @ epochs = {self.eq_epoch}
    AltAz @ epochs = {self.az_epoch}
"""
        return str
class ObservingSite(object):
    """Create an observing site

    Attributes:
        
        location: dictionary:
            lon: float [deg]: geodetic longitude
            lat: float [deg]: geodetic latitude
            elevation: float [km]: eleveation
        
        body: dictionary:
            id: SPICE id of body
            frameid: Rotating reference frame
    """
    def __init__(self,
                 location=dict(lon=0*u.deg,lat=0*u.deg,height=0),
                 conditions=dict(pressure=0*u.hPa,temperature=0*u.deg_C,relative_humidity=0,obswl=1*u.micron),
                 onplanet=None):

        # Body
        self.planet = onplanet
        
        # Geodetic coordinates
        self.location = location  

        # Atmospheric conditions
        self.conditions = conditions 

        # Cartesian coordinates
        #kwargs = location | conditions
        kwargs = location | dict()
        self.site = EarthLocation(**kwargs)
        self.ITRF = np.array(list(self.site.value))

        # Update site with a dummy time (J2000.0)
        self.update_site(MonTime(0,format='et',scale='tt'))

        # Update position with respect to body center
        self.geopos=spy.georec(self.location['lon'].value*DEG,
                               self.location['lat'].value*DEG,
                               self.location['height'].value*DEG,
                               self.planet.Re,self.planet.f)
        
        # From local refrence frame to body reference frame
        normal=spy.surfnm(self.planet.Re,
                          self.planet.Rt,
                          self.planet.Rp,
                          self.geopos)
        uy=spy.ucrss(np.array([0,0,1]),normal)
        ux=spy.ucrss(normal,uy)
        self.local2body=np.array(np.vstack((ux,uy,normal)).transpose().tolist())
        self.body2local=spy.invert(self.local2body)

        # Update orientation
        self.update_orientation()

    def update_site(self,mtime):  
        self.mtime = mtime
        et = mtime.et
        self.M_ITRF_J2000 = spy.pxform('IAU_EARTH','J2000',et)
        self.J2000 = spy.mxv(self.M_ITRF_J2000,np.array(list(self.site.value)))
        self.frame = AltAz(obstime=mtime.astrotime,location=self.site)
        self.update_orientation(et,framesky='EARTHTRUEEPOCH')

    def update_orientation(self,et=0,framesky='J2000'):
        """Update orientation matrices for an observing site

        Parameters:
            et: float [ephemeris seconds]:
                Time.

            framesky: string: default = 'J2000':
                Target with respect to which the transformation is performed
                Examples: J2000, ECLIPJ2000, EARTHTRUEEPOCH

        Update:
            Matrices:
                body2celestial, celestial2body.
        """
        #Local to Body reference frame transformation
        self.body2celestial=spy.pxform(self.planet.framebody,framesky,et)
        self.celestial2body=spy.invert(self.body2celestial)

    def to_celestial(self,az,el):
        """Convert from Azimuth and elevation to celestial coordinates
        """
        rlocal=spy.latrec(1,az*DEG,el*DEG)
        rbody=spy.mxv(self.local2body,rlocal)
        rcelestial=spy.mxv(self.body2celestial,rbody)
        r,RA,Dec=spy.reclat(rcelestial)
        RA=RA+2*np.pi if RA<0 else RA
        return RA*RAD,Dec*RAD
    
    def to_local(self,RA,Dec):
        """Convert from celestial coordinates to local coordinates
        """
        rcelestial=spy.latrec(1,RA*DEG,Dec*DEG)
        rbody=spy.mxv(self.celestial2body,rcelestial)
        rlocal=spy.mxv(self.body2local,rbody)
        r,az,el=spy.reclat(rlocal)
        az=az+2*np.pi if az<0 else az
        return az*RAD,el*RAD


    def __str__(self):
        str = f"""
Geodetic coordinates: {self.location}
Atmospheric conditions: {self.conditions}
Cartesian coordinates:
  ITRF = {self.ITRF}
  J2000 = {self.J2000}   
Planet: {self.planet}
Present epoch: {self.mtime.datestr}
Frame: {self.frame}
"""
        return str