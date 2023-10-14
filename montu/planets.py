from montu import *

###############################################################
# Planetary body
###############################################################
class PlanetaryBody_x(object):
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
        # exec(f'self.pyplanet = pyplanets_{self.capital}')

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
            Montu.vprint(verbose,"\tCoordinates @ J2000: ")
            Montu.vprint(verbose,"\t\tEquatorial:",Montu.dec2hex(self.RAJ2000),Montu.dec2hex(self.DecJ2000))
            Montu.vprint(verbose,"\t\tEcliptic:",Montu.dec2hex(self.LonJ2000),Montu.dec2hex(self.LatJ2000))
            Montu.vprint(verbose,f"\tCoordinates @ Epoch : ")
            Montu.vprint(verbose,"\t\tEquatorial:",Montu.dec2hex(self.RAEpoch),Montu.dec2hex(self.DecEpoch))
            Montu.vprint(verbose,"\t\tEcliptic:",Montu.dec2hex(self.LonEpoch),Montu.dec2hex(self.LatEpoch))
            Montu.vprint(verbose,f"\tLocal true sidereal time: ",Montu.dec2hex(self.tsa))
            Montu.vprint(verbose,f"\tHour angle @ Epoch: ",Montu.dec2hex(self.HA))
            Montu.vprint(verbose,f"\tLocal coordinates @ Epoch: ",Montu.dec2hex(self.az),Montu.dec2hex(self.el))

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
            Montu.vprint(verbose,"\tCoordinates @ J2000: ")
            Montu.vprint(verbose,"\t\tEquatorial:",Montu.dec2hex(self.RAJ2000),Montu.dec2hex(self.DecJ2000))
            Montu.vprint(verbose,"\t\tEcliptic:",Montu.dec2hex(self.LonJ2000),Montu.dec2hex(self.LatJ2000))
            Montu.vprint(verbose,f"\tCoordinates @ Epoch : ")
            Montu.vprint(verbose,"\t\tEquatorial:",Montu.dec2hex(self.RAEpoch),Montu.dec2hex(self.DecEpoch))
            Montu.vprint(verbose,"\t\tEcliptic:",Montu.dec2hex(self.LonEpoch),Montu.dec2hex(self.LatEpoch))
            Montu.vprint(verbose,f"\tLocal true sidereal time: ",Montu.dec2hex(self.tsa))
            Montu.vprint(verbose,f"\tHour angle @ Epoch: ",Montu.dec2hex(self.HA))
            Montu.vprint(verbose,f"\tLocal coordinates @ Epoch: ",Montu.dec2hex(self.az),Montu.dec2hex(self.el))

            compute = True

        if method == 'SPICE' or all_methods:

            # Using SPICE+pyplanets tools
            Montu.vprint(verbose,"Method 'SPICE':")

            # Retrieve positions in space
            site_planet_SSB_J2000,lt = spy.spkezr(site.planet.id,mtime.tt,'J2000','None','SSB')
            planet_SSB_J2000,lt = spy.spkezr(self.id,mtime.tt,'J2000','None','SSB')
            site_SSB_J2000 = site_planet_SSB_J2000[:3] + site.pos_J2000 

            # Celestial Coordinates at J2000
            planet_site_J2000 = planet_SSB_J2000[:3] - site_SSB_J2000
            r,RAJ2000,DECJ2000 = spy.recrad(planet_site_J2000)
            self.RAJ2000 = RAJ2000*RAD/15
            self.DecJ2000 = DECJ2000*RAD

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
            Montu.vprint(verbose,"\tCoordinates @ J2000: ")
            Montu.vprint(verbose,"\t\tEquatorial:",Montu.dec2hex(self.RAJ2000),Montu.dec2hex(self.DecJ2000))
            Montu.vprint(verbose,"\t\tEcliptic:",Montu.dec2hex(self.LonJ2000),Montu.dec2hex(self.LatJ2000))
            Montu.vprint(verbose,f"\tCoordinates @ Epoch : ")
            Montu.vprint(verbose,"\t\tEquatorial:",Montu.dec2hex(self.RAEpoch),Montu.dec2hex(self.DecEpoch))
            Montu.vprint(verbose,"\t\tEcliptic:",Montu.dec2hex(self.LonEpoch),Montu.dec2hex(self.LatEpoch))
            Montu.vprint(verbose,f"\tLocal true sidereal time: ",Montu.dec2hex(site.ltst))
            Montu.vprint(verbose,f"\tHour angle @ Epoch: ",Montu.dec2hex(self.HA))
            Montu.vprint(verbose,f"\tLocal coordinates @ Epoch: ",Montu.dec2hex(self.az),Montu.dec2hex(self.el))

            compute = True

        if not compute:
            raise ValueError(f"Method '{method}' for computing ephemerides not recognized")

        if store:
            self.df=pd.concat([self.df,pd.DataFrame(self.data)])
            self.df.reset_index(drop=True,inplace=True)

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
                # Method
                f'RAJ2000_{metstr}':[self.RAJ2000],f'DecJ2000_{metstr}':[self.DecJ2000],
                f'RAEpoch_{metstr}':[self.RAEpoch],f'DecEpoch_{metstr}':[self.DecEpoch],
                f'LonJ2000_{metstr}':[self.LonJ2000],f'LatJ2000_{metstr}':[self.LatJ2000],
                f'LonEpoch_{metstr}':[self.LonEpoch],f'LatEpoch_{metstr}':[self.LatEpoch],
                f'tsa_{metstr}':[self.tsa],f'HA_{metstr}':[self.HA],f'az_{metstr}':[self.az],f'el_{metstr}':[self.el],
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
class ObservingSite_x(object):
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
