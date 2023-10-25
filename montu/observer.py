from montu import *
###############################################################
# Module constants
###############################################################

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
        self.ltst = self.epoch.gtst + self.lon/15
